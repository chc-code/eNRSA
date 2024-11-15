#! /usr/bin/env python3

"""
functions used for NRSA package
"""
import os
import sys
import re
import time
import json, pickle
import gzip
from subprocess import Popen, PIPE
import pandas as pd
import numpy as np
import fisher
import traceback
import bisect 
import gc
from types import SimpleNamespace
from statsmodels.stats.multitest import multipletests
import statsmodels.stats.contingency_tables as contingency_tables
from scipy.stats import chi2
import warnings
import urllib.request
import hashlib
import tempfile
sys.dont_write_bytecode = True

pw_code = os.path.dirname(os.path.realpath(__file__))
inf_neg = float('-inf')
bases_set = set('ATGCatgc')
warnings.filterwarnings("ignore", category=RuntimeWarning, message="divide by zero encountered in scalar divide")

time_cost_util = {
    'bedtools_around_enhancer_center': 0,
    'bedtools_around_tss_for_heatmap': 0,
    
    }

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        self._original_stderr = sys.stderr
        sys.stdout = open(os.devnull, 'w')
        sys.stderr = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stderr.close()
        sys.stdout = self._original_stdout
        sys.stderr = self._original_stderr
def red(s, p=False):
    s_new = f'\u001b[31;1m{s}\u001b[0m'
    if p:
        print(s_new)
    return s_new
def green(s, p=False):
    s_new = f'\u001b[32;1m{s}\u001b[0m'
    if p:
        print(s_new)
    return s_new

def force_symlink(src, dst):
    try:
        os.symlink(src, dst)
    except FileExistsError:
        # If the symlink already exists, remove it and create a new one
        os.remove(dst)
        os.symlink(src, dst)

def getlogger(fn_log=None, logger_name=None, nocolor=False, verbose=False):
    import logging
    logger_name = logger_name or "main"
    
    try:
        logger = logging.getLogger(logger_name)
    except:
        logger = logging.getLogger('terminal')

    class CustomFormatter(logging.Formatter):
    
        def __init__(self, nocolor=False):
            self.nocolor = nocolor
        colors = {
            'black': '\u001b[30;1m',
            'red': '\u001b[31;1m',
            'r': '\u001b[31;1m',
            'bold_red': '\u001b[31;1m',
            'rb': '\u001b[31;1m',
            'green': '\u001b[32;1m',
            'g': '\u001b[32;1m',
            'gb': '\u001b[32;1m',
            'yellow': '\u001b[33;1m',
            'blue': '\u001b[34;1m',
            'b': '\u001b[34;1m',
            'purple': '\u001b[35;1m',
            'p': '\u001b[35;1m',
            'grey': '\u001b[38;1m',
        }
        FORMATS = {
            logging.WARNING: colors['purple'],
            logging.ERROR: colors['bold_red'],
            logging.CRITICAL: colors['bold_red'],
        }
    
        def format(self, record):
            format_str = "%(asctime)s  %(levelname)-6s %(funcName)-20s  line: %(lineno)-5s  %(message)s"
            reset = "\u001b[0m"
            log_fmt = None
            
            record.msg = str(record.msg)
            if self.nocolor:
                pass
            elif '@' in record.msg[:10]:
                try:
                    icolor, tmp = record.msg.split('@', 1)
                    log_fmt = self.colors.get(icolor)
                    if log_fmt:
                        record.msg = tmp
                except:
                    raise
                    pass
            else:
                log_fmt = self.FORMATS.get(record.levelno)
            if log_fmt:
                record.msg = log_fmt + record.msg + reset
            formatter = logging.Formatter(format_str, datefmt='%Y-%m-%d %H:%M:%S')
            return formatter.format(record)
    
    logger.setLevel('DEBUG')
    handler_names = {_.name for _ in logger.handlers}
    if 'console' not in handler_names:
        console = logging.StreamHandler(sys.stdout)
        console.setFormatter(CustomFormatter(nocolor=nocolor))
        console.setLevel('DEBUG' if verbose else 'INFO')
        console.name = 'console'
        logger.addHandler(console)

    if fn_log and 'file' not in handler_names:
        fh_file = logging.FileHandler(fn_log, mode='w', encoding='utf8')
        fh_file.setLevel('DEBUG')
        fh_file.setFormatter(CustomFormatter())
        fh_file.name = 'file'
        logger.addHandler(fh_file)
    return logger

logger = getlogger(logger_name='NRSA')

if os.path.exists('/.dockerenv'):
    in_docker = True
else:
    in_docker = False


def refine_chr(chr_):
    chr_ = chr_.strip().lower()
    for _ in ['chromosome', 'chrom', 'chr']:
        chr_ = chr_.replace(_, '')
    return {'23': 'x', '24': 'y', '25': 'mt', 'm': 'mt'}.get(chr_) or chr_

def build_idx_for_fa(home, fn_fa):
    # >chr1
    # NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    # NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    # NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    res = {'note': 'the line_width does not include the end newline symbol, it is the count of the regular char only'}
    line_width = {}
    # line_width,  the most frequent sequence in each line. the width does not include the end newline symbol
    # other keys = each chr, removed the chr prefix,  the value is like [s_byte, length]
    if '.gz' in fn_fa:
        logger.error(f'Please use decompressed fasta file')
        sys.exit(1)
        return None

    lb = os.path.basename(fn_fa).rsplit('.', 1)[0]
    fn_fa = os.path.abspath(fn_fa)
    fa_size = os.path.getsize(fn_fa)
    fa_name_hash = get_md5_str(fn_fa)[:8]
    
    
    # check if have write permission to the folder of the fa file
    if not os.access(os.path.dirname(fn_fa), os.W_OK):
        fa_pwout = f'{home}/.nrsa'
        os.makedirs(fa_pwout, exist_ok=True)
    else:
        fa_pwout = os.path.dirname(fn_fa)
        
    fn_idx = f'{fa_pwout}/{lb}.size_{fa_size}_{fa_name_hash}.idx.pkl'
    logger.debug(f'fa_pwout = {fa_pwout}, fn_idx = {fn_idx}')
    if os.path.exists(fn_idx):
        logger.debug(f'loading fasta index file: {fn_idx}')
        with open(fn_idx, 'rb') as f:
            res = pickle.load(f)
            return res

    logger.info(f'building index for {fn_fa}')
    chr_ = None
    with open(fn_fa) as f:
        while True:
            i = f.readline()
            if not i:
                break
            if i[0] == '>':
                chr_ = re.split(r'[\s\|:\[\]]+', i[1:].strip(), 1)[0] # avoid the chromosome with trailing strings
                chr_ = refine_chr(chr_)
                s_pos = f.tell()
                res[chr_] = [s_pos, 0]
            elif chr_:
                len1 = len(i) - 1
                res[chr_][1] += len1
                line_width.setdefault(len1, 0)
                line_width[len1] += 1
    
    tmp = sorted(line_width.items(), key=lambda _: _[1], reverse=True)
    word_wrap, freq = tmp[0]
    logger.debug(f'idx saved to {fn_idx}, word wrap = {word_wrap}, lines width = {tmp}')
    res['line_width'] = word_wrap
    with open(fn_idx, 'wb') as o:
        pickle.dump(res, o)
    return res

def dummy(idx, fh, chr_, start, end):
    return 'ATCG'

def get_seqence_from_fa(idx, fh, chr_, start, end):
    """
    idx is from loading of the json {fn_fa}.idx.json file, fh is the opened file handle of the fasta file
    already converted to uppercase
    for 20k transcripts, get the whole gene sequence cost about 2 second
    """
    # get_seqence_from_fa.note = 'normal'
    
    chr_ = refine_chr(chr_)
    line_width = idx['line_width']
    len1 = end - start + 1

    if chr_ not in idx:
        logger.error(f'invalid chromosome {chr_}, valid = {sorted(idx.keys())}')
        return 1
    byte_0, chr_total_len = idx[chr_]
    if end > chr_total_len:
        logger.error(f'end position out of range! chr{chr_} : total length = {chr_total_len} ')
        return 1
    lines_start = start // line_width
    lines_end = end // line_width
    lines_span = lines_end - lines_start
    start_byte = start + byte_0 -1 + lines_start
    fh.seek(start_byte)
    return fh.read(len1 + lines_span).replace('\n', '').upper()


def check_is_sorted(fn_bed):
    """
    check if the bed file is sorted or not, 
    """
    chr_prev = None
    chr_set = set()
    logger.info(f'Checking if bed file is sorted: {fn_bed}')
    chr_order = []
    with gzip.open(fn_bed, 'rt') if fn_bed.endswith('.gz') else open(fn_bed) as f:
        last = 0
        for i in f:
            chr_, start = i.split('\t', 2)[:2]
            start = int(start)
            if chr_ != chr_prev:
                chr_order.append(chr_)
                if chr_ in chr_set:
                    # the chr_ has been seen before
                    return False
                chr_set.add(chr_)
                chr_prev = chr_
            elif start < last:
                return False
            last = start
    # check if the chr is ordered
    tmp = sorted(chr_order)
    if tmp == chr_order:
        return True
    return False
def get_ref(home, organism, fn_gtf=None, fn_fa=None):
    # get the path of the reference files
    ref_files = {}
    pw_code = os.path.dirname(os.path.realpath(__file__))
    # check if the folder is writable
    if not os.access(pw_code, os.W_OK) or pw_code.startswith('/app/nrsa'):
        logger.debug(f'no write permission to the folder of the code: {pw_code}, will write to {home}/.nrsa')
        pw_code_write = f'{home}/.nrsa'
        os.makedirs(pw_code, exist_ok=True)
    else:
        pw_code_write = pw_code
    
    logger.debug(f'pw_code = {pw_code}, pw_code_write = {pw_code_write}, pw_code writeable = {os.access(pw_code, os.W_OK)}')
    if fn_gtf:
        fn_gtf = os.path.abspath(fn_gtf)
    if fn_fa:
        fn_fa = os.path.abspath(fn_fa)

    # create the folders if not exists
    for i in ['ref', 'fa', 'annotation']:
        os.makedirs(f'{pw_code_write}/{i}', exist_ok=True)
        
    # if the organism is not in the list, check fn_gtf and fn_fa, if either is None, return None, otherwise, use it and create a folder in pw_code/ref and pw_code/fa and create symlink
    builtin_organisms = ['ce10', 'dm3', 'dm6', 'hg19', 'hg38', 'mm10', 'mm39', 'danrer10']
    fn_extra_org = f'{pw_code_write}/ref/extra_organisms.txt'
    
    if os.path.exists(fn_extra_org):
        with open(fn_extra_org) as f:
            organism_extra = [_.strip().lower() for _ in f if _.strip()]
            builtin_organisms += organism_extra
    
    organism_lower = organism.lower()
    fn_gtf_exp = f'{pw_code}/ref/{organism}/RefSeq-{organism}.gtf'
    fn_fa_exp = f'{pw_code}/fa/{organism}/{organism}.fa'
    fn_gtf_exp_write = f'{pw_code_write}/ref/{organism}/RefSeq-{organism}.gtf'
    fn_fa_exp_write = f'{pw_code_write}/fa/{organism}/{organism}.fa'
    
    fn_protein_coding = f'{pw_code}/ref/{organism}/{organism}.protein_coding.txt'
    fn_protein_coding_write = f'{pw_code_write}/ref/{organism}/{organism}.protein_coding.txt'
    
    if fn_fa and not os.path.exists(fn_fa_exp_write):
        os.makedirs(os.path.dirname(fn_fa_exp_write), exist_ok=True)
        try:
            force_symlink(fn_fa, fn_fa_exp_write)
        except:
            pass
    if fn_gtf and not os.path.exists(fn_gtf_exp_write):
        os.makedirs(os.path.dirname(fn_gtf_exp_write), exist_ok=True)
        try:
            force_symlink(fn_gtf, fn_gtf_exp_write)
        except:
            pass
    
    ref_files = {}
    need_download = {}
    if organism_lower not in builtin_organisms:
        if fn_gtf is None or fn_fa is None:
            logger.error(f'New organism: {organism}, please provide the gtf and fa file by -gtf and -fa options')
            sys.exit(1)
        err = 0
        if not os.path.exists(fn_gtf):
            logger.error(f'GTF file not exist: {fn_gtf}')
            err = 1
        if not os.path.exists(fn_fa):
            logger.error(f'FA file not exist: {fn_fa}')
            err = 1
        if err:
            sys.exit(1)
        organism = organism_lower
        # create the folders if not exists and create the symlink
        fn_gtf = os.path.realpath(fn_gtf)
        fn_fa = os.path.realpath(fn_fa)
        
        logger.warning(f'new organism input, if you want to reuse this gtf / fa combination by specify -organism option, you can add the organism name as a new line in {fn_extra_org} and create symlink of gtf and fa files to\n{fn_gtf_exp_write}\n{fn_fa_exp_write}')

        ref_files['gtf'] = fn_gtf
        ref_files['fa'] = fn_fa
        ref_files['protein_coding'] = None
    else:
        # get the gtf and fa file of known genome
        for k, flist in [['gtf', [fn_gtf, fn_gtf_exp, fn_gtf_exp_write]], ['fa', [fn_fa, fn_fa_exp, fn_fa_exp_write]], ['protein_coding', [fn_protein_coding, fn_protein_coding_write]]]:
            found = 0
            for fn in flist:
                if fn is None:
                    continue
                if os.path.exists(fn):
                    ref_files[k] = fn
                    found = 1
                    break
            if not found:
                tmp = [(fn, os.path.exists(str(fn))) for fn in flist]
                logger.debug(f'{k} not found: {flist} - {tmp}')
                need_download[k] = flist[1]
        
    ref_files["fantom"], ref_files['association'] = {
        'hg19': ["human_permissive_enhancers_phase_1_and_2.bed", "human_enhancer_tss_associations.bed"],
        'mm10': ['mouse_permissive_enhancers_phase_1_and_2.bed', None],
        'hg38': ['human_permissive_enhancers_phase_1_and_2-hg38.bed', 'human_enhancer_tss_associations-hg38.bed']
        
    }.get(organism, [None, None])
    
    # 4d, dm3, dm6, hg19, hg38, mm10
    # chr1	557489	560146	chr5	134284878	134293544	PCBD2;	MCF7	ChIA-PET	19890323
    if organism in ['dm3', 'dm6', 'hg19', 'hg38', 'mm10']:
        ref_files['4d'] = f'{pw_code}/annotation/4DGenome-{organism}.txt'
    else:
        ref_files['4d'] = None
    
    if ref_files['fantom']:
        ref_files["fantom"] = f'{pw_code}/annotation/{ref_files["fantom"]}'
        if ref_files['association']:
            ref_files['association'] = f'{pw_code}/annotation/{ref_files["association"]}'

    # validite file existence
    for key in ref_files:
        value = ref_files[key]
        if not value:
            continue
        if not os.path.exists(value):
            if os.path.exists(f'{value}.gz'):
                logger.debug(f'use gzip version: {value}.gz')
                ref_files[key] = f'{value}.gz'
                continue
            # if exist in pw_code_write
            fn_new = value.replace(f'{pw_code}/', f'{pw_code_write}/')
            if pw_code != pw_code_write and os.path.exists(fn_new):
                logger.debug(f'exist in {pw_code_write}, but not in default ref folder {pw_code}: {fn_new}')
                ref_files[key] = fn_new
                continue
            tmp = value.split('/')
            if len(tmp) > 1 and tmp[-2] == 'annotation':
                if organism in {'mm10', 'hg19', 'hg38'}:
                    logger.warning(f"annotation file not exist: {organism} - {key}: '{value}'\n")
                    need_download[key] = value
                ref_files[key] = None
            else:
                logger.warning(f"{organism} - {key}: '{value}'  not exist\n")
                need_download[key] = value
    need_download = {k: v.replace(f'{pw_code}/', f'{pw_code_write}/') for k, v in need_download.items()}
    if need_download:
        logger.warning(f'{len(need_download)} reference files need to be downloaded for {organism}...')
        logger.debug(need_download)
        folders = {os.path.dirname(fn) for fn in need_download.values()}
        logger.debug(f'folders needed = {folders}')
        for ipw in folders:
            os.makedirs(ipw, exist_ok=True)
        ref_files = download_ref_files(home, organism, need_download, ref_files)
    else:
        logger.debug(f'All reference files exist for {organism}')

    tmp = json.dumps(ref_files, indent=3)
    logger.debug(tmp)
    return ref_files

def get_md5(fn):
    with open(fn, "rb") as f:
        md5 = hashlib.md5()
        while chunk := f.read(8192):
            md5.update(chunk)
    return md5.hexdigest()

def get_md5_str(s):
    md5 = hashlib.md5()
    md5.update(s.encode())
    return md5.hexdigest()

def download_file_task(url, save_path, md5_exp=None):
    try:
        urllib.request.urlretrieve(url, save_path)
        # check the md5sum
        if md5_exp:
            logger.info(f'Calculating md5 for {save_path}')
            md5_actual = get_md5(save_path)
            if md5_actual != md5_exp:
                logger.error(f"MD5 not match: {save_path}, expected = {md5_exp}, actual = {md5_actual}")
                return 1
        file_size = os.path.getsize(save_path)
        logger.info(f"File downloaded successfully and saved to {save_path}, size = {file_size/1024/1024}MB")
        return 0
    except:
        e = traceback.format_exc()
        logger.error(f"Error downloading the file - {save_path}, url = {url}: \n{e}")
        return 1

def download_ref_files(home, organism, need_download, ref_files):
    """
    download the missing reference files
    """
    builtin_organisms = ['ce10', 'dm3', 'dm6', 'hg19', 'hg38', 'mm10', 'mm39', 'danrer10']
    if organism not in builtin_organisms:
        logger.error(f'organism {organism} not in the builtin list, please provide the gtf and fa file')
        sys.exit(1)
    url_base = 'https://bioinfo.vanderbilt.edu/NRSA/download/NRSA'
    err = 0
    
    url_md5sum = f'{url_base}/files_md5sum.json'
    if not os.access(pw_code, os.W_OK) or pw_code.startswith('/app/nrsa'):
        logger.debug(f'no write permission to the folder of the code: {pw_code}, will write to {home}/.nrsa')
        pw_code_write = f'{home}/.nrsa'
        os.makedirs(pw_code, exist_ok=True)
        os.makedirs(pw_code + '/ref', exist_ok=True)
    else:
        pw_code_write = pw_code    
    
    fn_md5sum_local = f'{pw_code_write}/ref/files_md5sum.json'
    status = download_file_task(url_md5sum, fn_md5sum_local)
    if status:
        logger.warning(f'MD5 checksum file not downloaded: {url_md5sum}')
        md5sum_all = {}
    else:
        with open(fn_md5sum_local) as f:
            md5sum_all = json.load(f)
    
    for key, fn_dest in need_download.items():
        slug = {'gtf': organism, 'fa': 'fa', 'fantom': 'annotation', 'association': 'annotation', '4d': 'annotation'}[key]
        base_name = os.path.basename(fn_dest)
        md5_exp = md5sum_all.get(base_name)
        url_full = f'{url_base}/{slug}/{base_name}'
        status = download_file_task(url_full, fn_dest, md5_exp)
        if status:
            err = 1
        else:
            ref_files[key] = fn_dest
    if err:
        sys.exit(1)
    return ref_files

def bench(s, lb, d):
    now = time.time()
    d[lb] += now - s
    return now


def get_peak_method1(count_per_base, chr_, strand_idx, s, e):
    # sum from the per-base pre-count dict
    ct_sum = 0
    # strand_idx = 0 if strand == '+' else 1
    return sum([count_per_base.get(i, [0, 0])[strand_idx] for i in range(s, e + 1)])

def get_peak_method2(count_per_base, count_bin, chr_, strand_idx, s, e, bin_size):
    # compared to get_peak_method1 function
    # test a 1.5kb region to get the mapped reads count
    # method2, combination of per-site count and bin_collected count (bin_size = 10), cost 19us 
    # method1, use the per-site count only, cost 150us
    # so the method2 is 7 times faster than method1
    # for small window, e.g 50bp, method2 is still faster, 5.8us vs 10.8us
    # will also test if will improve when the pre_count bin_size is 20bp
    # for bin_size = 20, speed is 12.5 and 5.4us for 1.5kb and 50bp window, both faster than bin_size = 10
    
    
    # res is outer dict, stores the per_site count
    # bin_size = count_bin['bin_size']
    bin_start = s // bin_size
    bin_end = e // bin_size
    if bin_end - bin_start < 2:
        return sum([count_per_base.get(i, [0, 0])[strand_idx] for i in range(s, e + 1)])
    # the region is covered by 3 or more bins
    points = []
    left_mod = s % bin_size
    right_mod = e % bin_size
    if left_mod:
        points += range(s, (bin_start + 1) * bin_size)
    if right_mod != bin_size - 1:
        bin_list = range(bin_start + 1, bin_end)
        points += range(bin_end * bin_size, e + 1)
    else:
        bin_list = range(bin_start + 1, bin_end + 1)
    return sum([count_per_base.get(i, [0, 0])[strand_idx] for i in points]) + sum([count_bin.get(i, [0, 0])[strand_idx] for i in bin_list])

def get_peak(count_per_base, count_bin, chr_, strand, gene_raw_s, strand_idx, pp_start, pp_end, gb_start, gb_end, tts_start, tts_end, gb_len_mappable, gene_seq, window_size, step_size, bin_size, prev_pp_peak):
    # {'chr': '5', 'strand': '-', 'gene_name': 'PFDN1', 'start': 139682626, 'end': 139682689}
    # pro_up = pro_down = 500
    # gb_start_pos = 1000 # count the gb density from 1000bp downstream of the gene start site
    # using the already_parsed_windows will cost more time

    gbc = get_peak_method2(count_per_base, count_bin, chr_, strand_idx, gb_start, gb_end, bin_size=bin_size)
    if gbc == 0:
        gbd = 0
    else:
        gbd = gbc / gb_len_mappable
    # gbc = gbd = 0
    # k_pp_region = f'{chr_}{strand}:{pp_start}'
    
    # get tts_count
    tts_ct = get_peak_method2(count_per_base, count_bin, chr_, strand_idx, tts_start, tts_end, bin_size=bin_size)
    
    
    # get all the sites count in pp region
    pp_region_count = [count_per_base.get(i, [0, 0])[strand_idx] for i in range(pp_start, pp_end + 1)]
    window_ct = [(i, sum(pp_region_count[i: i + window_size])) for i in range(0, pp_end - pp_start - window_size + 2, step_size)]
    peak_window_start, ppc  = max(window_ct, key=lambda _: _[1]) 
    peak_window_ct = [(i, pp_region_count[i]) for i in range(peak_window_start, peak_window_start + window_size)]
    summit_pos, summit_count = max(peak_window_ct, key=lambda _: _[1])
    summit_pos += pp_start
    summit_pos_str = f'{chr_}{strand}:{summit_pos}'
    
    if ppc == 0:
        ppd = 0
    elif gene_seq:
        peak_windows_start_abs = peak_window_start + pp_start - gene_raw_s
        mappable_sites = window_size - gene_seq[peak_windows_start_abs: peak_windows_start_abs + window_size].count('N')
        ppd = ppc / mappable_sites
    else:
        ppd = ppc / window_size

    pp_res = {'ppc': ppc, 'ppd': ppd, 'mappable_sites': window_size, 'summit_pos': summit_pos_str, 'summit_count': summit_count}
    # prev_pp_peak[k_pp_region] = pp_res
    # gbc = gbd = 0
    # return 0, 0, {'ppc': 0, 'ppd': 0, 'mappable_sites': 50, 'summit_pos': 0, 'summit_count': 0}
    return gbc, gbd,  pp_res, tts_ct


def draw_box_plot(n_gene_cols, pwout, out_name, n_rep1, n_rep2=None, gn_list=None):
    """
    gn_list, the genes to be plotted, if None, all genes will be plotted
    pwout: folder name for pause_PROseq.pl outputs
    n_rep1 and n_rep2:  the replacates number for condition1 and conditon2
    n_rep2 is optional
    """
    from matplotlib import pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    
    n_rep2 = n_rep2 or 0
    gn_list = set(gn_list) if gn_list else None
    # check the pause_PROseq.pl output files
    for i in ['known_gene', 'intermediate']:
        ipw = f'{pwout}/{i}'
        if not os.path.exists(ipw):
            logger.error(f'output folder of pause_PROseq.py not found: {i}')
            return 1

    # read pp_gb results
    fn_pp_gb = f'{pwout}/known_gene/normalized_pp_gb.txt'
    fn_pindex = f'{pwout}/known_gene/pindex.txt'
    fno = f'{pwout}/intermediate/boxplot_data.txt'
    
    header_str = []
    out_str = {} # key = gn, v = list combined from pp_gb and pindex
    
    n_total_sam = n_rep1 + n_rep2
    
    # parse pp_gb results
    with open(fn_pp_gb) as f:
        header_in = f.readline().strip().split('\t')
        if header_in[0] != 'Transcript':
            logger.error(f'invalid header in {fn_pp_gb}')
            return 1
        idx_gene_cols = list(range(n_gene_cols))
        idx_list = []
        idx_list += idx_gene_cols
        idx_list += [i for i, col in enumerate(header_in) if col[:4] == 'ppd_']
        idx_list += [i for i, col in enumerate(header_in) if col[:4] == 'gbd_']
        header_str = [header_in[i] for i in idx_list]

        for i in f:
            line = i.strip().split('\t')
            gn = line[0]
            if gn_list is None or gn in gn_list:
                k = f'{gn}@{line[1]}'
                ires = [line[i] for i in idx_list]
                out_str[k] = ires

    # parse pindex results
    with open(fn_pindex) as f:
        header_in = f.readline().strip().split('\t')
        # transcript, gene, sam1-pindex, sam1-pvalue, sam1-FDR, sam2-pindex, sam2-pvalue, sam2-FDR...
        if header_in[0] != 'Transcript':
            logger.error(f'invalid header in {fn_pindex}')
            return 1
        idx_list = [i for i, col in enumerate(header_in) if col.endswith('pindex')]
        header_str += [header_in[_] for _ in idx_list]
        
        for i in f:
            line = i[:-1].split('\t')
            gn = line[0]
            if gn_list is None or gn in gn_list:
                k = f'{gn}@{line[1]}'
                if k not in out_str:
                    logger.warning(f'{gn} in pindex not found in pp_gb')
                    continue
                ires = out_str[k]
                try:
                    vals_pindex = [line[_] for _ in idx_list]
                except:
                    logger.error(line)
                    sys.exit(1)
                if '' in vals_pindex:
                    logger.error(line)
                    sys.exit(1)
                ires += vals_pindex
                

    # write to file
    logger.debug(f'data for box plot saved to {fno}')
    with open(fno, 'w') as o:
        o.write('\t'.join(header_str) + '\n')
        for k in sorted(out_str):
            o.write('\t'.join(out_str[k]) + '\n')
    
    # draw the box plot, adapted from boxplot.R
    # R CMD --args outdir=\"$inter_dir\" pname=\"$name\" custom=1 rep1=$rep1 rep2=$rep2
    df_box = pd.read_csv(fno, sep='\t')
    cols_pp_density = [_ for _ in df_box.columns if _.startswith('ppd_')]
    cols_gb_density = [_ for _ in df_box.columns if _.startswith('gbd_')]
    cols_pindex = [_ for _ in df_box.columns if _.endswith('pindex')]
    sam_list = [name.split('ppd_')[1] for name in cols_pp_density]
    # logger.info(df_box.head())
    
    pp_density = df_box.loc[:, cols_pp_density]
    gb_density = df_box.loc[:, cols_gb_density]
    pindex =     df_box.loc[:, cols_pindex]
    
    # logger.info(pindex.head(20))

    pp_density.columns = sam_list
    gb_density.columns = sam_list
    pindex.columns = sam_list
    
    # plot
    fn_plot_pp = f'{pwout}/known_gene/{out_name}_pp_density.pdf'
    fn_plot_gb = f'{pwout}/known_gene/{out_name}_gb_density.pdf'
    fn_plot_pindex = f'{pwout}/known_gene/{out_name}_pausing_index.pdf'
    
    def plot_task(fn, ylabel, data, factor=1):
        with PdfPages(fn) as pdf:
            plt.figure(figsize=(2, 4))
            
            # Split the data into two groups based on n_rep1
            data_rep1 = data.iloc[:, :n_rep1]
            data_rep1 = [data_rep1[_].dropna() * factor for _ in data_rep1]
            
            if n_rep2:
                data_rep2 = data.iloc[:, n_rep1:]
                data_rep2 = [data_rep2[_].dropna() * factor for _ in data_rep2]
            
            # previously, data_rep1 and data_rep2 were all x 1000 as the input of box plot data
            
            # Plot the boxplot for the first n_rep1 samples (blue color)
            # plt.boxplot(data_rep1 * 1000, positions=range(1, n_rep1 + 1), showfliers=False, patch_artist=True, boxprops=dict(facecolor='blue'), capprops=dict(color='blue'), whiskerprops=dict(color='blue'), medianprops=dict(color='blue'))
            plt.boxplot(data_rep1, positions=range(1, n_rep1 + 1), showfliers=False, patch_artist=True, boxprops=dict(facecolor='blue'), capprops=dict(color='blue'), whiskerprops=dict(color='blue'), medianprops=dict(color='blue'))
            
            # Plot the boxplot for the rest of the samples (red color)
            if n_rep2:
                # plt.boxplot(data_rep2 * 1000, positions=range(n_rep1 + 1, n_total_sam + 1), showfliers=False, patch_artist=True, boxprops=dict(facecolor='red'), capprops=dict(color='red'), whiskerprops=dict(color='red'), medianprops=dict(color='red'))
                plt.boxplot(data_rep2, positions=range(n_rep1 + 1, n_total_sam + 1), showfliers=False, patch_artist=True, boxprops=dict(facecolor='red'), capprops=dict(color='red'), whiskerprops=dict(color='red'), medianprops=dict(color='red'))
            
            plt.ylabel(ylabel)
            plt.xticks(range(1, n_total_sam + 1), sam_list, rotation='vertical')
            plt.tight_layout()
            pdf.savefig()
            
    plot_task(fn_plot_pp, 'Reads per Kb (RPK)', pp_density, 1000)
    plot_task(fn_plot_gb, 'Reads per Kb (RPK)', gb_density, 1000)
    plot_task(fn_plot_pindex, 'Pausing Index', pindex)

    def plot_hist(n_sam, n_prev_col, ppd, gbd, condition_sn):

        plt.figure(figsize=(12, 6 * (n_sam - 1)))
        for k in range(1, n_sam):
            col1 = n_prev_col
            col2 = k + n_prev_col
            col1_str = ppd.columns[col1]
            col2_str = ppd.columns[col2]
            xlabel = f"log2({sam_list[0]}/{sam_list[k]})"
            plt.subplot((n_sam - 1), 2, k)
            ppd_tmp = ppd.loc[(ppd.iloc[:, col1] > 0) & (ppd.iloc[:, col2] > 0), [col1_str, col2_str]]
            plt.hist(np.log2(ppd_tmp.iloc[:, 0] / ppd_tmp.iloc[:, 1]), bins=100)
            plt.title(f"ppd: rep1vs.rep{k + 1}")
            plt.xlabel(xlabel)

            plt.subplot((n_sam - 1), 2, k + 1)
            gbd_tmp = gbd.loc[(gbd.iloc[:, col1] > 0) & (gbd.iloc[:, col2] > 0), [col1_str, col2_str]]
            plt.hist(np.log2(gbd_tmp.iloc[:, 0] / gbd_tmp.iloc[:, 1]), bins=100)

            plt.title(f"gbd: rep1vs.rep{k + 1}")
            plt.xlabel(xlabel)

        plt.savefig(f"{pwout}/known_gene/Reps-condition{condition_sn}.tif")
        plt.close()
    
    if n_rep1 > 1:
        plot_hist(n_rep1, 0, pp_density, gb_density, 1)
    if n_rep2 > 1:
        plot_hist(n_rep2, n_rep1, pp_density, gb_density, 2)
    

def draw_heatmap_pindex(pwout, fdr_thres=0.05, fc_thres=0):
    """
    plot for the pindex change
    fc_thres, the fc absolution value threshold for the pindex change
    """
    from matplotlib import pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import LinearSegmentedColormap
    
    fn_pindex_change = f'{pwout}/known_gene/pindex_change.txt'
    fno = f'{pwout}/known_gene/pindex_change.pdf'
    if not os.path.exists(fn_pindex_change):
        logger.error(f'{fn_pindex_change} not found')
        return 1

    data_plot = pd.read_csv(fn_pindex_change, sep='\t')
    # Transcript	Gene	log2fc	pvalue	FDR
    col_fdr, col_logfc = 'FDR', 'log2fc'
    data_plot[col_logfc] = data_plot[col_logfc].replace([np.inf, -np.inf], np.nan)

    # drop the rows which have NA in the log2fc and fdr column
    data_plot = data_plot.dropna(subset=[col_logfc, col_fdr])

    # filter out the genes with FDR > 0.05
    data_plot = data_plot[(data_plot[col_fdr] < fdr_thres) & (abs(data_plot[col_logfc]) > fc_thres)]
    # sort by log2FC
    data_plot = data_plot.sort_values(by=col_logfc)

    with PdfPages(fno) as pdf:
        plt.figure(figsize=(2, 6))
        grid = plt.GridSpec(2, 1, height_ratios=[14, 1])
        ax1 = plt.subplot(grid[0, 0])
        values = data_plot['log2fc'].values[::-1]
        cmap = LinearSegmentedColormap.from_list("my_palette", ["green", "yellow", "red"], N=209)
        ax1.imshow(values.reshape(-1, 1), cmap=cmap, aspect='auto')
        ax1.axis('off')
        
        
        ax2 = plt.subplot(grid[1, 0])
        gradient = np.linspace(min(values), max(values), 209).reshape(1, -1)
        ax2.imshow(gradient, aspect='auto', cmap=cmap)
        ax2.set_xticks([0, 104, 208])
        ax2.set_xticklabels([round(min(values), 2), "0", round(max(values), 2)], fontsize=8)
        ax2.set_yticks([])
        ax2.set_xlabel("")
        ax2.set_ylabel("")
        plt.subplots_adjust(left = 0, right = 0.99, top=0.95, bottom = 0, hspace=0.05)
        pdf.savefig(bbox_inches='tight')
        plt.close()
    
def get_line_count(fn):
    n = 0
    with gzip.open(fn) if fn.endswith('.gz') else open(fn) as f:
        for i in f:
            n += 1
    return n

def calculate_signal(fn_lb, fn_bed, fn_split_bin_bed, fn_coverage, demo=False):
    if demo and os.path.exists(fn_coverage):
        logger.debug(f'skipped bedtools coverage for enhancer center, {fn_lb}')
    else:
        logger.debug(f'bedtools coverage around enhancer center: {fn_lb}')
        s = time.time()
        cmd = f'bedtools coverage -a {fn_split_bin_bed} -b {fn_bed} -counts -sorted > {fn_coverage}'
        retcode, stdout, stderr = run_shell(cmd, ret_info=True)
        time_cost_util['bedtools_around_enhancer_center'] += time.time() - s
        if retcode:
            logger.error(f'failed to run bedtools coverage for {fn_lb}')
            logger.debug(stderr)
            return 1
    # process the coverage file
    # chr1	323932	323952	1   10
    count = {}
    with open(fn_coverage) as f:
        for i in f:
            line = i[:-1].split('\t')
            bin_sn, ict = line[3], line[4]
            count.setdefault(bin_sn, 0)
            count[bin_sn] += int(ict)
    return {int(k): v for k, v in count.items()}

def draw_signal(pwout, fn_enhancer_center, fls_ctrl, fls_case, chr_map, distance=2000, bin_size=20, signal_type='p', demo=False):
    # fn_enhancer_center #enhancer central position file, the output of eRNA.pl (Enhancer_centralposition.txt); required, enhancer central position file, the output of eRNA.pl (Enhancer_center.txt), or the file of chromation location of interest (should contains at least two columns separated by tab or space, e.g "chr11	83078317")
    # example line = chr1	565340	565340	+, no header
    
    # my $distance = 2000;	#up and down relative distance from peak centre;
    # my $bin_size=20; #bin size for smoothing;
    # my @case_bed;	#read alignment files in bed/bam format (case);
    # my @control_bed;	#read alignment files in bed/bam format (control);
    # my $signal_type="p"; The signal_type (-s) should be either 'p' or 'o', p - signal for PROseq was input in pause_PROseq.pl; o - signal for other data such as histone modification (default: p)

    # fn_enhancer_center is like below, with header
    # chr	position	FONTOM5	Enhancer_ID
    # 1	8215322	N	95
    # 1	8299373	N	96
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    fls_case = fls_case or []
    if signal_type not in {'p', 'o'}:
        logger.error(f'invalid signal type: {signal_type}, valid = "p" or "o"')
        return 1
    pwerna = f'{pwout}/eRNA'
    pwinter = f'{pwout}/intermediate'
    fn_padding = f'{pwinter}/enhancer_center.padding.tmp'
    fn_split_bin_bed = f'{pwinter}/enhancer_center.split_regions.bed'

    n_valid_sites = 0
    n_bad_line = 0 # the line which 2nd column is not numeric
    min_dist_before = distance + bin_size/2
    nfactor = {}  # only used when signal_type= "o",  key = fn_lb, v = total lines in the bed file, the line_count should be calculated during the pre-counting step
    res = {}
    half_bin_size = int(bin_size / 2)
    

    with open(fn_enhancer_center) as f, open(fn_padding, 'w') as o:
        f.readline()  # skip header
        err = None
        for i in f:
            try:
                line = i[:-1].split('\t')
                pos = int(line[1])
                if pos > min_dist_before:
                    n_valid_sites += 1
                    chr_ = chr_map.get(line[0], line[0])
                    left_boundary = pos - distance - half_bin_size
                    right_boundary = pos + distance + half_bin_size
                    print(f'{chr_}\t{left_boundary}\t{right_boundary}', file=o)
            except:
                if err is None:
                    err = traceback.format_exc()
                n_bad_line += 1
    if n_bad_line:
        logger.warning(f'{n_bad_line} lines in {fn_enhancer_center} are invalid, 2nd column is not numeric')
        logger.debug(f'the error info is {err}')

    # get the line count for each bed file
    for fn_lb, fn_bed in fls_ctrl + fls_case:
        fn_line_count = fn_bed.replace('.sorted.bed', '.line_count.txt')
        if os.path.exists(fn_line_count):
            with open(fn_line_count) as f:
                nfactor[fn_lb] = int(f.readline())
        elif signal_type == 'o':
            line_count = get_line_count(fn)
            nfactor[fn_lb] = line_count
            with open(fn_line_count, 'w') as o:
                o.write(str(line_count))

    # build the regions around the enhancer center for plot
    # no region name, the final output from makewindows is like
    # chr1	20	30	3  the 4th column is the region number
    def makewindows(fn_lb, fn_tss_padding, bin_size, chr_map):
        # rewrite the bedtools makewindows, because the original one does not contain the strand information
        fno = f'{pwout}/intermediate/{fn_lb}.split_regions.bed'
        with open(fn_tss_padding) as f, open(fno, 'w') as o:
            for i in f:
                chr_, s, e = i[:-1].split('\t')
                bin_sn = 0
                chr_orig = chr_map.get(chr_, chr_)
                s, e = int(s), int(e)
                for pos in range(s, e, bin_size):
                    bin_sn += 1
                    print(f'{chr_orig}\t{pos}\t{pos + bin_size}\t{bin_sn}', file=o)
        # sort it
        retcode = run_shell(f'bedtools sort -i {fno} > {fno}.tmp && mv {fno}.tmp {fno}')
        return fno
    
    # build the new tss for the first file, and use it  for the rest
    fn_split_bin_bed = makewindows(fn_lb, fn_padding, bin_size, chr_map)

    # get the signal results
    for fn_lb, fn_bed in fls_ctrl + fls_case:
        fn_coverage = f'{pwinter}/{fn_lb}.signal.coverage.txt'
        isignal = calculate_signal(fn_lb, fn_bed, fn_split_bin_bed, fn_coverage, demo=demo)
        if isignal == 1:
            logger.error(f'Fail to process signal for {fn_lb} around the enhancer centers')
            return 1
        res[fn_lb] = isignal
    
    # save nfactor for type = 'o'
    if signal_type == 'o':
        with open(f'{pwinter}/nf_o.txt') as o:
            print('sample\tnfactor', file=o)
            for fn_lb, line_count in nfactor.items():
                print(f'{fn_lb}\t{10_000_000/line_count}', file=o)

    # export the signal results
    fn_pro_signal = f'{pwinter}/PROseq_signal.txt'
    sam_list = [fn_lb for fn_lb, _ in fls_ctrl + fls_case]
    
    with open(fn_pro_signal, 'w') as o:
        print('position\t' + '\t'.join(sam_list), file=o)
        bin_sn = 0
        for bin_pos in range(-distance, distance + 1, bin_size):
            bin_sn += 1
            row = [str(bin_pos)]
            for fn_lb, _ in fls_ctrl + fls_case:
                row.append(f'{res[fn_lb][bin_sn]/n_valid_sites:.4f}')
            print('\t'.join(row), file=o)
    


    logger.debug('plot signal')
    fn_signal_pdf = f'{pwout}/eRNA/signal_around_ehancer-center.pdf'

    fn_nf = f'{pwout}/intermediate/' + (f'nf.txt' if signal_type == 'p' else 'nf_o.txt')
    def process_nf_line(line):
        lb, v = line.strip().split('\t')
        return lb, float(v)
    with open(fn_nf) as f:
        f.readline() # skip header
        nf = dict([process_nf_line(i) for i in f if i.strip()])
    nf_array = np.array([nf[fn_lb] for fn_lb in sam_list])
    
    
    data = pd.read_csv(fn_pro_signal, sep='\t')
    # get the normalized count
    data[sam_list] = data[sam_list] * nf_array

    n_ctrl, n_case = len(fls_ctrl), len(fls_case)

    color_map = {
        'ctrl': LinearSegmentedColormap.from_list("blue_darkblue", ["blue", "darkblue"], N=n_ctrl),
        'case': LinearSegmentedColormap.from_list("red_darkred", ["red", "darkred"], N=n_case)
    }
    
    plt.figure(figsize=(6, 6))
    linestyle_pool = ['-', ':', '--', '-.']
    len_line_styles = len(linestyle_pool)
    for condition, fls in zip(['ctrl', 'case'], [fls_ctrl, fls_case]):
        colors = color_map[condition]
        sam_idx = 0
        n_sam_condition = len(fls)
        for fn_lb, _ in fls:
            if n_sam_condition > 1:
                normalized_index = sam_idx / (n_sam_condition - 1)
            else:
                normalized_index = 0
            linestyle = linestyle_pool[sam_idx % len_line_styles]
            plt.plot(data['position'], data[fn_lb], color=colors(normalized_index), linestyle=linestyle, label=fn_lb, linewidth=0.9)
            sam_idx += 1
            
    # Adding legend
    plt.legend(sam_list, loc='upper right', fontsize=9, bbox_to_anchor=(1, 1), borderaxespad=0.)

    # Labels and plot configuration
    plt.xlabel('Distance to enhancer center', fontsize=15)
    plt.ylabel('Fragment per bp per peak', fontsize=15)
    plt.ylim(0, data.iloc[:, 1:].max().max() * 1.2)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)

    # Save plot as PDF
    plt.savefig(fn_signal_pdf, format='pdf', bbox_inches='tight')
    plt.close()


    # plot the signal
    # signal = 1 if signal_type == 'o' else 0
    # file=fn_pro_signal  outdir=\"$inter_dir\" sample=\"$samplename_str\" pname=\"$name\" rep1=$rep1 rep2=$rep2 signal=1 \' $cmd_file

def draw_heatmap_pp_change(n_gene_cols, pwout, pw_bed, fls_ctrl, fls_case, fn_tss, region_size=5000, bin_size=200, outname='heatmap', skipe_bedtools_coverage=False):
    """
    draw heatmap for pp change
    fls_ctrl, fls_case, foramt is like [fn_lb, fn_bed]
    region_size, upstream and downstream distance relative to TSS for plotting PRO-seq signal (bp, default: 5000),should can be divided by bin size
    bin_size (bp, default: 200),should can be divided by 2
    """
    from matplotlib import pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.backends.backend_pdf import PdfPages

    fls_case = fls_case or []
    
    ppchange = {}
    fn_pp_change = f'{pwout}/known_gene/pp_change.txt'
    fn_out_pdf = f'{pwout}/known_gene/heatmap.pdf'
    # fn_data_bed = f'{pwout}/intermediate/data_bed.tmp'
    fn_tss_padding = f'{pwout}/intermediate/infile_bed.tmp' # expand the fn_tss with upstream and downstream region by region_size
    time_bedtools_coverage = 0

    if bin_size % 2 or region_size % bin_size:
        logger.error(f'bin_size should be even and region_size should be divisible by bin_size')
        return 1

    bin_number = int(region_size * 2 / bin_size + 1)  # *2 because of the upstream and downstream

    def to_number(s):
        if not s or s.upper() == 'NA':
            return 'NA'
        if 'inf' in s:
            return 'NA'
        return float(s)

    # read pp_change results
    # Transcript	Gene	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
    n_case, n_ctrl = len(fls_case), len(fls_ctrl)
    
    # gene_list is the genes with non-NA log2FoldChange in the pp_change.txt
    with open(fn_pp_change) as f:
        header = f.readline().strip().split('\t')
        logger.debug(header)
        idx_transcript = header.index('Transcript')
        idx_logfc = header.index('log2FoldChange')
        for i in f:
            line = i[:-1].split('\t')
            transcript_id = line[idx_transcript]
            logfc = to_number(line[idx_logfc])
            if logfc != 'NA':
                ppchange[transcript_id] = logfc

    # tss file
    # chr1	11874	11874	NR_046018.2	+
    # tss_res = []
    tss_with_padding_n = 0
    strand_info = {}
    half_bin_size = int(bin_size / 2)
    with open(fn_tss) as f, open(fn_tss_padding, 'w') as o:
        for i in f:
            line = i.strip().split('\t')
            chr_, pos, transcript_id, strand = [line[_] for _ in [0, 1, 3, 4]]
            pos_int = int(pos)
            chr_ = refine_chr(chr_)
            if transcript_id in ppchange and ppchange[transcript_id] != 'NA':
                logfc = ppchange[transcript_id]
                row = [transcript_id, chr_, pos, strand, logfc]
                strand_info[transcript_id] = strand

                left_boundary = pos_int - region_size - half_bin_size
                right_boundary = pos_int + region_size + half_bin_size
                if left_boundary < 0:
                    continue
                row_tmp = [chr_, left_boundary, right_boundary, transcript_id, strand]
                print('\t'.join(map(str, row_tmp)), file=o)
                # tss_res.append(row)
                tss_with_padding_n += 1
    if tss_with_padding_n == 0:
        logger.error(f'no gene in tss file match with the gene list, please check if the transcript names in the GTF file match with the TSS file')
        return 1
    
    def makewindows(fn_lb, fn_tss_padding, bin_size, chr_map):
        # rewrite the bedtools makewindows, because the original one does not contain the strand information
        fno = f'{pwout}/intermediate/{fn_lb}.split_regions.bed'
        with open(fn_tss_padding) as f, open(fno, 'w') as o:
            for i in f:
                chr_, s, e, ts, strand = i[:-1].split('\t')
                s, e = int(s), int(e)
                bin_sn = 0
                chr_orig = chr_map.get(chr_, chr_)
                for pos in range(s, e, bin_size):
                    bin_sn += 1
                    print(f'{chr_orig}\t{pos}\t{pos + bin_size}\t{ts}_{bin_sn}\t-1\t{strand}', file=o)
        # logger.info('sort split region bed 2')
        fntmp = f'{fno}.tmp'
        retcode = run_shell(f'bedtools sort -i {fno} > {fntmp} && mv {fntmp} {fno}')

        return fno
    
    # build the new tss for the first file, and use it  for the rest
    fn_lb = fls_ctrl[0][0]
    fn_chr_map = f'{pw_bed}/{fn_lb}.chr_map.pkl'
    with open(fn_chr_map, 'rb') as f:
        chr_map = pickle.load(f)
    split_bed_per_file = makewindows(fn_lb, fn_tss_padding, bin_size, chr_map)
    
    # load the normalization factors (nf.txt)
    fn_nf = f'{pwout}/intermediate/nf.txt'
    factors_dict = {}
    with open(fn_nf) as f:
        f.readline() # header
        for i in f:
            sample, factor = i.strip().split('\t')
            factors_dict[sample] = float(factor)
    
    # logger.debug(factors_dict)
    # main part, get the overlap of the input bed files with the split bin bed
    # coverage_by_strand_flag = '-s'  # by_strand, if used, need to update the reference tss file
    # logger.warning(f'currently getting coverage by strand')

    coverage_by_strand_flag = '' # current setting
    for condition, fls in {'control': fls_ctrl, 'case': fls_case}.items():
        fn_count_sum = f'{pwout}/intermediate/{condition}.count'
        count = {}   # the sum of norm count of each bin across all samples in this condition
        n_sam_condition = len(fls)
        for fn_lb, fn_bed  in fls:
            norm_factor = factors_dict[fn_lb]
            
            fn_coverage_tmp = f'{pwout}/intermediate/{fn_lb}.coverage_count.tmp'
            # do we need to calculate the coverage by strand??? 

            # if not :
            if skipe_bedtools_coverage and os.path.exists(fn_coverage_tmp):
                logger.debug(f'skip bedtools coverage for {fn_lb}')
            else:
                logger.info(f'getting coverage for {fn_lb}')
                s = time.time()
                cmd = f'bedtools coverage -a {split_bed_per_file} -b {fn_bed} -counts {coverage_by_strand_flag} -sorted > {fn_coverage_tmp}'
                logger.debug(cmd)
                retcode, msg_stdout, msg_stderr = run_shell(cmd, ret_info=True)
                if retcode:
                    logger.error(f'failed to run bedtools coverage for {fn_lb}')
                    logger.debug(msg_stderr)
                    return 1
                # process the coverage file
                # chr1	323932	323952	NR_028322.1_8   10
                dure = time.time() - s
                time_cost_util['bedtools_around_tss_for_heatmap'] += dure
                time_bedtools_coverage += dure
                logger.debug(f'done, time used = {dure:.2f}s')


            with open(fn_coverage_tmp) as f:
                for i in f:
                    try:
                        # chr13	51846175	51846375	NM_011817.2_24	+	9
                        _, transcript_id_chunk, _, strand, ict = i.strip().rsplit('\t', 4)
                        if ict == '0':
                            continue
                        transcript_id, bin_sn = transcript_id_chunk.rsplit('_', 1)
                        # bin_sn = int(bin_sn)
                        ict = int(ict)
                        ict_norm = ict * norm_factor
                        count.setdefault(transcript_id, {}).setdefault(bin_sn, 0)
                        count[transcript_id][bin_sn] += ict_norm
                    except:
                        e = traceback.format_exc()
                        logger.error(e)
                        logger.error(f'error when processing {fn_coverage_tmp}, line = \n{i}\n')
                        sys.exit(1)

        with open(fn_count_sum, 'w') as o:
            half_bin_number = (bin_number - 1) // 2
            # header = ['Transcript'] + [f'up_{i}' for i in range(half_bin_number, 0, -1)] + ['tss'] + [f'down_{i}' for i in range(1, half_bin_number + 1)]
            header = [f'up_{i}' for i in range(half_bin_number, 0, -1)] + ['tss'] + [f'down_{i}' for i in range(1, half_bin_number + 1)]
            print('\t'.join(header), file=o)
            bin_num_plus = [str(_) for _ in range(1, bin_number + 1)]
            bin_num_minus = [str(_) for _ in range(bin_number, 0, -1)]
            for transcript_id, v1 in count.items():
                strand = strand_info[transcript_id]
                bin_count_list = []
                bin_number_order = bin_num_plus if strand == '+' else bin_num_minus
                
                for bin_sn in bin_number_order:
                    bin_count_list.append(str(round(v1.get(bin_sn, 0) / n_sam_condition, 2))) # average the count across all samples
                tmp = "\t".join(bin_count_list)
                print(f'{transcript_id}\t{tmp}', file=o)
    
    logger.debug(f'total time used for bedtools coverage: {time_bedtools_coverage:.2f}s')

    # heatmap.R  --args file=\"$data_bed\" outdir=\"$inter_dir\" pname=\"$outname\" window=$bin_size region=$region_size
    # each row is [transcript_id, chr_, pos, strand, logfc]
    # read the count table
    df_case = pd.read_csv(f'{pwout}/intermediate/case.count', sep='\t')
    df_ctrl = pd.read_csv(f'{pwout}/intermediate/control.count', sep='\t')
    case_log = np.log2(df_case + 1)
    ctrl_log = np.log2(df_ctrl + 1)
    df_delta = case_log - ctrl_log
    abs_values = np.abs(df_delta.values.flatten())
    cutoff = np.quantile(abs_values, 0.75) # verified, match with R code
    cutoff1 = round(cutoff, 1)
    df_delta = df_delta.clip(lower=-cutoff, upper=cutoff)
    # logger.info(df_delta.head())
    
    
    # plot
    n_color_block = 5
    with PdfPages(fn_out_pdf) as pdf:
    # Set up the layout and figure size
        plt.figure(figsize=(6, 18))
        grid = plt.GridSpec(2, 1, height_ratios=[12, 1])

        lab = round(region_size / 1000, 1) # label for the heatmap
        # Define color palette
        colors = ["blue", "yellow"]

        # Create a colormap with 5 discrete colors
        my_palette = LinearSegmentedColormap.from_list("custom_palette", colors, N=n_color_block)
        n_cols = df_delta.shape[1]
        
        ax1 = plt.subplot(grid[0, 0])
        ax1.imshow(df_delta, cmap=my_palette, vmin=-cutoff, vmax=cutoff, aspect='auto')
        ax1.set_xticks([0, n_cols // 2, n_cols - 1])
        ax1.set_xticklabels([f"-{lab}K", "0", f"+{lab}K"], fontsize=12)
        ax1.set_yticks([])
        ax1.set_xlabel("")
        ax1.set_ylabel("")
        
        
        ax2 = plt.subplot(grid[1, 0])
        ax2.imshow(np.arange(1, n_color_block + 1).reshape(1, -1), cmap=my_palette, aspect='auto')
        ax2.set_xticks([0, 2, 4])
        ax2.set_xticklabels([f'-{cutoff1}', '0', f'+{cutoff1}'])
        ax2.set_yticks([])
        ax2.set_xlabel("")
        ax2.set_ylabel("")
        
        
        pdf.savefig(bbox_inches='tight')
        plt.close()

def groseq_fisherexact_pausing_for_gene(c1, l1, c2, l2):
    """
    Perform Fisher's exact test to calculate the p-value for the pausing index of each gene. x is a row of the pausing count dataset
    per run, is about 1ms, so it is quite slow., 1-2ms per run
    user new fisher package, the speed is 100x faster (6us per call)
    must make sure the input are all non-zero
    the pvalue is slightly different from R fisher.test(x, alternative='two.sided') esp. for small pvalues
    """
    if not all([c1, l1, c2, l2]):
        return np.nan
    # null hypothesis: the pause read count is uniformly distributed in the gene body and the promoter region
    expectC1 = round(l1 * (c1 + c2) / (l1 + l2))
    expectC2 = c1 + c2 - expectC1
    # default = two-sided
    return fisher.pvalue(c1, expectC1, c2, expectC2).two_tail

def groseq_fisherexact_comparison_for_gene(tssCountGene1, gbCountGene1, tssCountGene2, gbCountGene2):
    """
    Perform Fisher's exact test to calculate the p-value for the comparison of pausing index between two genes or 2 conditions of the same gene
    """
    c1 = round(float(tssCountGene1))
    c2 = round(float(gbCountGene1))
    c3 = round(float(tssCountGene2))
    c4 = round(float(gbCountGene2))
    
    if c1 == 0 or c2 == 0 or c3 == 0 or c4 == 0:
        return np.nan
    
    # default = two-sided
    return fisher.pvalue(c1, c2, c3, c4).two_tail


def get_FDR_per_sample(fn_lb, fn_peak, fn_fdr, pause_index_str):
    """
    get the pvalue and FDR based on the peak file, retain only 1 transcript per gene, (keep the transcript with the highest pp count)
    """
    logger.debug('loading raw peak table')
    # df_peak = pd.read_csv(fn_peak, sep='\t')
    df_peak_gene = pd.read_csv(fn_peak, sep='\t')
    
    # below code is for collapse genes
    # but actually, do not collapse them. this will lead to lots of problems for combining FDR files from different samples. e.g. building the pindex.txt
    # Transcript	Gene	ppc	ppm	ppd	pps	gbc	gbm	gbd	pauseIndex
    # NM_000014.5	A2M	40	50	0.8	12-:0	43	47522	0.0009	884.13023
    # idx = df_peak.groupby('Gene')['ppc'].idxmax()
    # df_peak_gene = df_peak.loc[idx]
    logger.debug('get pvalue by fisher exact test')
    df_peak_gene['pvalue'] = df_peak_gene.apply(lambda x: groseq_fisherexact_pausing_for_gene(x['ppc'], x['ppm'], x['gbc'], x['gbm']), axis=1)
    
    pvals = df_peak_gene['pvalue'].dropna()
    logger.debug('get FDR')
    fdr = multipletests(pvals, method='fdr_bh')[1]
    df_peak_gene.loc[pvals.index, 'FDR'] = fdr
    logger.debug(f'dump to file {fn_fdr}')
    df_peak_gene.to_csv(fn_fdr, sep='\t', index=False, na_rep='NA')
    
    logger.debug('building pause_index_str')
    col_ts = df_peak_gene.columns[0]
    df_tmp= df_peak_gene[[col_ts, 'pauseIndex', 'pvalue', 'FDR']].fillna('NA').astype(str)
    for row in df_tmp.itertuples(index=False, name=None):
        ts, *values = row
        pause_index_str.setdefault(ts, {})[fn_lb] = values
    return pause_index_str


def run_deseq2(n_gene_cols, data, metadata, ref_level, col_group=None, min_reads=10, filter_out_low_expression=False, 
               size_factors_in=None, size_factor_only=False, test_level=None):
    """
    A count matrix of shape ‘number of samples’ x ‘number of genes’, containing read counts (non-negative integers),
    Metadata (or “column” data) of shape ‘number of samples’ x ‘number of variables’, containing sample annotations that will be used to split the data in cohorts.  rows = sample, columns = gene
    index should be samplename
    meta data, first col is the sample name (index), following columns usually condidion and group
    min_reads:  if a gene has less than min_reads in all samples, it will be removed
    ref_level,  the level used as reference in the comparison
    """
    
    
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats
    from pydeseq2.preprocessing import deseq2_norm_fit
    
    # data filtering
    
    counttable = data.iloc[:, n_gene_cols:]
    idx_gene = 1
    idx_transcript = 0
    idx_gene_cols = list(range(n_gene_cols))
    # if n_gene_cols < 2:
    #     logger.error(f'DESeq analysis is for known genes only, the first 2 columns should be Transcript and Gene')

    # if each group only have one sample
    n_sam = counttable.shape[1]
    if n_sam == 2:
        log2fc = np.log2(data.iloc[:, n_gene_cols + 1]) / (data.iloc[:, n_gene_cols])
        res_df = data.iloc[:, idx_gene_cols]
        res_df['log2FoldChange'] = log2fc
        return res_df, np.array([1, 1])

    gene_col = data.iloc[:, idx_gene] # use Transcript as index
    counttable.index = data.iloc[:, idx_transcript]
    
    counttable = counttable.T
    metadata.index = counttable.index
    if 'batch' in metadata.columns:
        with_batch = True
    else:
        with_batch = False
    
    logger.debug(f'{list(metadata.columns)} , with_batch = {with_batch}')
    
    # remove samples for which annotations are missing and exclude genes with very low levels of expression
    col_condition = 'condition'
    design_factors = [col_condition, 'batch'] if with_batch else col_condition
    sample_to_keep = list(~metadata[col_condition].isna())
    counttable = counttable.loc[sample_to_keep, :]
    metadata = metadata.loc[sample_to_keep, :]
    
    if filter_out_low_expression:
        genes_to_keep = counttable.columns[counttable.sum(axis=0) >= min_reads]
        counttable = counttable[genes_to_keep]

    inference = DefaultInference(n_cpus=8)
    # logger.info(counttable.head())
    # logger.info(metadata)
    # logger.debug(data.shape)
    # logger.debug(data.head())
    # logger.debug(counttable.shape)
    # logger.debug(counttable.head())
    # logger.debug(metadata)
    # logger.debug(ref_level)
    levels = set(metadata[col_condition])
    levels.remove(ref_level)
    test_level_l = list(levels)
    
    if test_level:
        pass
    elif len(test_level_l) == 1:
        test_level = test_level_l[0]
    elif len(test_level_l) == 0:
        logger.error(f'no test level found in metadata, column = {col_condition}: \n{metadata}')
        return None, np.array([1, 1])
    elif len(test_level_l) > 0:
        logger.error(f'Test level not specified, column = {col_condition}: \n{metadata}')
        return None, np.array([1, 1])
    
    
    dds = DeseqDataSet(
        counts=counttable,
        metadata=metadata,
        design_factors=design_factors, # compare samples based on condition
        ref_level=[col_condition, ref_level],
        refit_cooks=True,
        inference=inference,
    )
    if size_factors_in is None:  # gb change
        # tested, the size factors are exactly the same as the R deseq2 code
        dds.fit_size_factors()
        size_factors = dds.obsm['size_factors']
    else:
        size_factors = size_factors_in
        dds.obsm["size_factors"] = size_factors

    if size_factor_only:
        return None, size_factors


    if size_factors_in is None:  # gb change
        dds.deseq2()
    else:
        # logger.debug(f'custom size factor specified: {size_factors_in}')
        counts = dds.X
        (dds.logmeans, dds.filtered_genes) = deseq2_norm_fit(counts)
        dds.layers["normed_counts"] = counts / size_factors[:, None]
        
        # added below for 0.4.10 version of pydeseq2
        dds.varm["_normed_means"] = dds.layers["normed_counts"].mean(0)


        dds.fit_genewise_dispersions()
        dds.fit_dispersion_trend()
        dds.fit_dispersion_prior()
        dds.fit_MAP_dispersions()
        dds.fit_LFC()
        dds.calculate_cooks()
        dds.refit()

    # avialble attributes = X: count tdata, obs = design factors, 
    # obsm = sample matrix, keys = design_matrix, size_factors
    # varm = gene matrix, keys = dispersions, LFC

    # get stat
    stat_res = DeseqStats(dds, inference=inference, contrast=[col_condition, test_level, ref_level])
    stat_res.summary()
    
    # the logFC is the same as R, but the p-value and padj are slightly different
    res_df = stat_res.results_df
    # here must convert gene_col to list, otherwise, the value will be all NaN, because the index does not match
    colname_transcript = data.columns[idx_transcript]
    res_df.reset_index(colname_transcript, inplace=True)
    
    # logger.debug('head of data')
    # logger.debug(data.head().to_dict())
    
    # logger.debug('header of res_df')
    # logger.debug(res_df.head().to_dict())
    
    res_df = data.iloc[:, :n_gene_cols].merge(res_df, on=colname_transcript, how='right')
    # res_df.insert(1, 'Gene', list(gene_col))
    res_df = res_df.sort_values('padj')
    
    return res_df, size_factors


def filter_pp_gb(data, n_gene_cols, rep1, rep2, skip_filtering=False):
    """
    reused by change_pp_gb_with_case and change_pindex
    """
    n_sam = rep1 + rep2
    cols_raw = data.columns
    col_idx = {
        'ppc': {
            'cond1': list(range(n_gene_cols, n_gene_cols + rep1)),
            'cond2': list(range(n_gene_cols + rep1, n_gene_cols + n_sam)),
            'combined': list(range(n_gene_cols, n_gene_cols + n_sam)),
        },
        'gbc': {
            'cond1': list(range(n_gene_cols + n_sam,             n_gene_cols + n_sam + rep1 * 2, 2)),
            'cond2': list(range(n_gene_cols + n_sam + rep1 * 2, n_gene_cols + n_sam * 3, 2)),
            'combined': list(range(n_gene_cols + n_sam, n_gene_cols + n_sam * 3, 2)),
        },
        'gbd': {
            'cond1': list(range(n_gene_cols + n_sam + 1,             n_gene_cols + n_sam + rep1 * 2, 2)),
            'cond2': list(range(n_gene_cols + n_sam + rep1 * 2 + 1,  n_gene_cols + n_sam * 3, 2)),
            'combined': list(range(n_gene_cols + n_sam + 1, n_gene_cols + n_sam * 3, 2)),
        },
    }
    
    # cond1 = ctrl, cond2 = case
    idx_ppc_combined = col_idx['ppc']['combined']
    idx_gbc_combined = col_idx['gbc']['combined']
    idx_gbd_combined = col_idx['gbd']['combined']
    idx_gene_cols = list(range(n_gene_cols))

    sam_list = [re.sub('^ppc_', '', _) for _ in data.columns[idx_ppc_combined]]
    
    if skip_filtering:
        data_pass = data
        data_drop = pd.DataFrame(columns=data.columns)
        return col_idx, sam_list, idx_ppc_combined, idx_gbc_combined, idx_gbd_combined, idx_gene_cols, data_pass, data_drop, None
    
    tmp1 = sum([1 if cols_raw[i].startswith('gbc_') else 0 for i in idx_gbc_combined])
    tmp2 = sum([1 if cols_raw[i].startswith('gbd_') else 0 for i in idx_gbd_combined])
    if tmp1 != n_sam or tmp2 != n_sam:
        logger.error('the number of samples in gbc or gbd is not correct')
        sys.exit(1)
        return None
    
    
    data['tmp_ppc_pass'] = data.iloc[:, idx_ppc_combined].apply(lambda x: all(x > 0), axis=1)
    
    sum_gbc = []
    for sn, idx_gbc_tmp, idx_gbd_tmp in zip(range(n_sam), idx_gbc_combined, idx_gbd_combined):
        sum_col = data.iloc[:, idx_gbc_tmp].sum()
        # logger.info(f'{sn} - {sum_col}')
        sum_gbc.append(int(sum_col))
        data[f'tmp_gbd_pass_{sn+1}'] = data.iloc[:, idx_gbd_tmp] * 10_000_000 / sum_col
    cols_gbd_pass = [f'tmp_gbd_pass_{sn+1}' for sn in range(n_sam)]
    data['tmp_gbd_pass'] = data[cols_gbd_pass].apply(lambda x: x.sum()/n_sam > 0.004, axis=1)

    
    # already compared with R result, the same
    data_pass = data[data['tmp_ppc_pass'] & data['tmp_gbd_pass']] # data0 in perl code
    data_drop = data.drop(data_pass.index) # data_0 in perl code

    data_pass = data_pass.drop(columns=[_ for _ in data_pass.columns if _[:4] == 'tmp_'])
    data_drop = data_drop.drop(columns=[_ for _ in data_drop.columns if _[:4] == 'tmp_'])

    return col_idx, sam_list, idx_ppc_combined, idx_gbc_combined, idx_gbd_combined, idx_gene_cols, data_pass, data_drop, sum_gbc

def change_pp_gb_with_case(n_gene_cols, rep1, rep2, data, pwout, window_size, factor_flag, factor1=None, factor2=None, islongerna=False):
    # data = processed, count_pp_gb.txt, already removed the chr, start, end  and strand column, and already collapsed the transcripts with the same gene
    # if factor1, factor2 is not None, it should be a list
    # n_gene_cols = 2 for known gene, and 1 for eRNA
    
    # columns = "Transcript\tGene"  + ppc_[sam_list], gbc_sam1, gbd_sam1, gbc_sam2, gbd_sam2 ....
    if rep2 == 0:
        factor2 = []
        
    if islongerna:
        pw_change_prefix = f'{pwout}/eRNA/longeRNA-'
    else:
        pw_change_prefix = f'{pwout}/known_gene/'

    n_sam = rep1 + rep2
    skip_filtering = islongerna # skip if islongerna
    col_idx, sam_list, idx_ppc_combined, idx_gbc_combined, idx_gbd_combined, idx_gene_cols, data_pass, data_drop, sum_gbc = filter_pp_gb(data, n_gene_cols, rep1, rep2, skip_filtering=skip_filtering)
    
    # logger.warning(n_gene_cols)
    # logger.warning(data.head())
    # logger.warning(data_pass.head())
    # logger.warning(data_drop.head())
    
    if not islongerna:
        data_pass.iloc[:, idx_gene_cols].to_csv(f'{pwout}/intermediate/active_gene.txt', sep='\t', index=False, header=False, na_rep='NA')
        
        fn_gbc_sum = f'{pwout}/intermediate/gbc_sum.json'
        with open(fn_gbc_sum, 'w')  as o:
            json.dump(sum_gbc, o)    
    data_pass_pp = data_pass.iloc[:, idx_gene_cols + idx_ppc_combined] # datapp in R
    data_drop_pp = data_drop.iloc[:, idx_gene_cols + idx_ppc_combined] # data_pp in R
    
    data_pass_gb = data_pass.iloc[:, idx_gene_cols + idx_gbc_combined] # datagb in R
    data_drop_gb = data_drop.iloc[:, idx_gene_cols + idx_gbc_combined] # data_gb in R
    
    fn_metadata = f'{pwout}/intermediate/design_table_for_deseq2.txt'
    metadata = pd.read_csv(fn_metadata, sep='\t', index_col=0)
    ref_level = 'control'
    test_level = 'case'


    # with normalization factors
    if factor_flag == 1:
        factor1 = np.array(factor1)
        factor2 = np.array(factor2)
        norm_factors = np.concatenate([factor1, factor2])
        size_factors = 1 / norm_factors
    else:  # without normalization factors
        size_factors = None
    
    # logger.warning(data_pass_gb.head())
    # logger.warning(data_drop_gb.head())
    # rows_with_na = data_pass_pp.loc[pd.isna(data_pass_pp).any(axis=1)]
    # logger.warning(rows_with_na.shape)
    # logger.warning(rows_with_na.head())
    
    if rep2 > 0:
        # gb change
        logger.info('running DESeq2')
        with HiddenPrints():
            res_df, size_factors = run_deseq2(n_gene_cols, data_pass_gb, metadata, ref_level, size_factors_in=size_factors, test_level=test_level)
            
            if res_df is None:
                sys.exit(1)
            
            res_df_full = pd.concat([res_df, data_drop_gb.iloc[:, idx_gene_cols]])
            
            res_df_full.to_csv(f'{pw_change_prefix}gb_change.txt', sep='\t', index=False, na_rep='NA')
            
            # pp change
            try:
                res_df, _ = run_deseq2(n_gene_cols, data_pass_pp, metadata, ref_level=ref_level, size_factors_in=size_factors, test_level=test_level)
            except:
                fn_deseq2_input = f'{pwout}/intermediate/deseq2.pkl'
                
                logger.error(f'fail to run deseq2, dumping the arguments to {fn_deseq2_input}')
                with open(fn_deseq2_input, 'wb') as o:
                    pickle.dump([n_gene_cols, data_pass_pp, metadata, ref_level, size_factors], o)
                
                sys.exit(1)
                
            res_df_full = pd.concat([res_df, data_drop_pp.iloc[:, idx_gene_cols]])
            res_df_full.to_csv(f'{pw_change_prefix}pp_change.txt', sep='\t', index=False, na_rep='NA')
    elif rep1 == 1:
        logger.debug('single contrl sample only')
        size_factors = 1 # size factor is 1 for the only sample
    elif rep2 == 0:
        # control only, multiple ctrol samples
        _, size_factors = run_deseq2(n_gene_cols, data_pass_gb, metadata, ref_level, size_factors_in=size_factors, size_factor_only=True, test_level=test_level)

    norm_factors = 1 / size_factors
    # get the normalized data in other function, because we need the chr start end strand information, and the data here is already processed
    # save normalization factors
    if not islongerna:
        nf = pd.DataFrame({'sample': sam_list, 'nfactor': norm_factors})
        nf.to_csv(f'{pwout}/intermediate/nf.txt', sep='\t', index=False, na_rep='NA')
    return size_factors, sam_list
    

def change_pp_gb(n_gene_cols, fn, pwout, rep1, rep2, window_size, factor1=None, factor2=None, factor_flag=0, islongerna=False):
    """
    adapted from Rscript change_pp_gb.R
    if factor_flag == 1, then will use these normalization factors, otherwise , will use DESeq2 to normalize the data
    fn = count_pp_gb.txt

    """
    # input file = my $out_pp_gb = $inter_dir . "count_pp_gb.txt";
    # columns = "Transcript\tGene\tchr\tstart\tend\tstrand"  + ppc_[sam_list], gbc_sam1, gbd_sam1, gbc_sam2, gbd_sam2 ....
    err = 0
    rep2 = rep2 or 0
    
    if factor1 is not None and not isinstance(factor1, list):
        logger.error('factor1 should be a list')
        err = 1
    if factor2 is not None and not isinstance(factor2, list):
        logger.error('factor2 should be a list')
        err = 1
    if err:
        return 1
    if islongerna:
        pw_change_prefix = f'{pwout}/eRNA/longeRNA-'
        n_extra_cols = 0 # chr, start, end, strand
    else:
        pw_change_prefix = f'{pwout}/known_gene/'
        n_extra_cols = 4 # chr, start, end, strand

    fn_norm = f'{pw_change_prefix}normalized_pp_gb.txt'

    data_raw = pd.read_csv(fn, sep='\t')
    data = data_raw.copy()
    data = data.iloc[:, list(range(n_gene_cols)) + list(range(n_gene_cols + n_extra_cols, len(data.columns)))]
    cols = data.columns
    cols_ppc = [i for i in cols if i.startswith('ppc_')]
    cols_gbc = [i for i in cols if i.startswith('gbc_')]
    cols_gbd = [i for i in cols if i.startswith('gbd_')]
    
    # only collapse the genes for refgenes, not longeRNA
    if not islongerna:
        data['ppc_sum'] = data[cols_ppc].sum(axis=1)
        data['gbc_sum'] = data[cols_gbc].sum(axis=1)
        
        ## !!! collapse the data from transcript level to gene level
        def is_max(col):
            max_v = col.max()
            return [1 if i == max_v else 0 for i in col]

        # colapase by ppc_sum
        data['SameGene_ppc_sum_max'] = data.groupby('Gene')['ppc_sum'].transform(is_max)
        data = data.loc[data['SameGene_ppc_sum_max'] == 1]
        
        # if multiple transcripts have the same max ppc_sum, then choose the one with the max gbc_sum
        data['SameGene_gbc_sum_max'] = data.groupby('Gene')['gbc_sum'].transform(is_max)
        data = data.loc[data['SameGene_gbc_sum_max'] == 1]
        
        # leaving the first one for each gene
        data = data.drop_duplicates(subset='Gene', keep='first').drop(columns=['SameGene_ppc_sum_max', 'SameGene_gbc_sum_max', 'ppc_sum', 'gbc_sum'])
    
    n_sam = rep1 + rep2
    
    # (rep1, rep2, data, pwout, window_size, factor_flag, factor1=None, factor2=None)
    size_factors, sam_list = change_pp_gb_with_case(n_gene_cols, rep1, rep2, data, pwout, window_size, factor_flag, factor1, factor2, islongerna=islongerna)

    # get the normalized data
    n_prev_cols = n_gene_cols + n_extra_cols
    col_idx, sam_list, idx_ppc_combined, idx_gbc_combined, idx_gbd_combined, idx_gene_cols, data_pass, data_drop, gbc_sum = filter_pp_gb(data_raw, n_prev_cols, rep1, rep2, skip_filtering=True)
    
    norm_factors = 1 / size_factors
    ppc_norm = data_raw.iloc[:, idx_ppc_combined] * norm_factors
    ppd_norm = ppc_norm / window_size
    ppd_norm.columns = [f'ppd_{_}' for _ in sam_list]
    
    gbc_norm = data_raw.iloc[:, idx_gbc_combined] * norm_factors
    gbd_norm = data_raw.iloc[:, idx_gbd_combined] * norm_factors
    data_normalized = pd.concat([data_raw.iloc[:, :n_prev_cols], ppc_norm, ppd_norm, gbc_norm, gbd_norm], axis=1)
    data_normalized.to_csv(fn_norm, sep='\t', index=False, na_rep='NA')


def cmhtest(arr):
    """
    performs the Mantel-Haenszel test on a certain gene
    the row is like ppc_sam1, ppc_sam2, ppc_sam3, ppc_sam4, gbc_sam1, gbd_sam1, gbc_sam2, gbd_sam2, ...
    """
    if len(arr) == 0:
        return 'NA', 'NA'
    arr_raw = arr
    arr = np.array(arr).T
    cmt = contingency_tables.StratifiedTable(tables=arr)
    odds_ratio = cmt.oddsratio_pooled
    test_res = cmt.test_null_odds(correction=True)
    stat_v = test_res.statistic
    pvalue = chi2.sf(stat_v, 1) # df = 1 for contigency table of 2x2
    # someof the pvalue is still 0, but it match with R, e.g.
    # [[[2894.0, 1340.0], [2504.0, 3180.0]], [[4566.0, 2048.0], [2504.0, 3180.0]], [[2894.0, 1340.0], [3262.0, 3275.0]], [[4566.0, 2048.0], [3262.0, 3275.0]]]
    # [2894.0, 1340.0, 2504.0, 3180.0, 4566.0, 2048.0, 2504.0, 3180.0, 2894.0, 1340.0, 3262.0, 3275.0, 4566.0, 2048.0, 3262.0, 3275.0]
    
    # statistic = test_res.statistic
    # print(f'odds_ratio = {odds_ratio}, pvalue = {pvalue}, statistic = {statistic}') 
    # return log2(odd_ratio), pvalue
    return odds_ratio if odds_ratio != 0 else 'NA', pvalue


def get_pvalue_2_sample(row, idx_ppc1, idx_ppc2, idx_gbc1, idx_gbc2):
    """
    get the p-value for 2 samples
    using groseq_fisherexact_comparison_for_gene(tssCountGene1, gbCountGene1, tssCountGene2, gbCountGene2)
    """
    tssCountGene1 = row[idx_ppc1]
    gbCountGene1 = row[idx_gbc1]
    tssCountGene2 = row[idx_ppc2]
    gbCountGene2 = row[idx_gbc2]
    pvalue = groseq_fisherexact_comparison_for_gene(tssCountGene1, gbCountGene1, tssCountGene2, gbCountGene2)
    if pvalue is None:
        return 'NA'
    return pvalue

def add_FDR_col(df, col_pvalue, col_FDR='FDR'):
    df = df.reset_index(drop=True)
    pd.set_option('future.no_silent_downcasting', True)
    pvals = df[col_pvalue].replace('NA', np.nan).astype(float).dropna()
    fdr = multipletests(pvals, method='fdr_bh')[1]
    # if 'Ctrl_max_count' in df.columns:
    #     logger.warning(f'df shape= {df.shape}, pvals shape = {pvals.shape}, fdr = {len(fdr)}')
    #     logger.debug(len(pvals.index))
    #     logger.debug(pvals.head(20))

    df.loc[pvals.index, col_FDR] = fdr
    fdr = df.pop(col_FDR)
    idx_pvalue = df.columns.get_loc(col_pvalue)
    df.insert(idx_pvalue + 1, col_FDR, fdr)
    return df


def change_pindex(fno_prefix, n_gene_cols, fn, fno, rep1, rep2, window_size, factor1=None, factor2=None, factor_flag=0):
    """
    fno_prefix, default is empty for known genes, longeRNA- for longeRNA
    adapted from Rscript change_pindex.R
    input = intermediate/ count_pp_gb.txt or longeRNA-count_pp_gb.txt  Transcript, Gene, chr, start, end, strand, (ppc_sam1 - ppc_sum2 ...), (gbc_sam1, gbd_sam1, gbc_sam2, gbd_sam2 ...)
    
    The transcripts were not collapsed to gene level
    """
    # for longeRNA, the count_pp_gb header is like
    # long-eRNA	ppc_ ...., so there is no Transcript, gene, chr, start, end, and strand
    # compared with previous results, matches
    if fno_prefix not in {'', 'longeRNA-'}:
        logger.error(f'invalid file prefix for pindex_change, should be empty or "longeRNA-", input="{fno_prefix}"')
        return 1
    if not rep2:
        logger.warning(f'Only control samples were specified, skip performing change_pindex...')
        return 
    is_erna = True if fno_prefix == 'longeRNA-' else False
    data = pd.read_csv(fn, sep='\t')
    n_sam = rep1 + rep2
    if is_erna:
        # the data does not contains the chr, start, end, strand columns, and there is no Gene column
        # no filtering needed
        ana_type = 'eRNA'
        n_cols_prev = n_gene_cols = 1
    else:
        ana_type = 'Known genes'
        n_cols_prev = n_gene_cols + 4
        cols_keep = list(range(n_gene_cols)) + list(range(n_cols_prev, len(data.columns))) # drop the chr, start, end, strand columns
        data = data.iloc[:, cols_keep]

    col_idx, sam_list, idx_ppc_combined, idx_gbc_combined, idx_gbd_combined, idx_gene_cols, data_pass, data_drop, gbc_sum = filter_pp_gb(data, n_gene_cols, rep1, rep2, skip_filtering=is_erna)

    def cmhtest_for_pindex_change(row):
        arr = []
        row = list(row)
        n_sam = rep1 + rep2
        for i in range(rep1):
            pp1 = row[i] # pp condition 1
            gb1 = row[n_sam + (i + 1) * 2 - 2] # gb condition 1
            for j in range(rep2):
                pp2 = row[rep1 + j]
                gb2 = row[n_sam + rep1 * 2 + (j + 1) * 2 - 2]
                # below is the gb2 defined in R code, which is wrong, here just used to test the consistency
                # gb2 = row[(rep2 + rep1 + j + 1) * 2 - 2]
                if pp1 + pp2 + gb1 + gb2 > 0 and gb2 * pp1 > 0:
                    arr.append([[pp2, gb2], [pp1, gb1]])
        odds_ratio, pvalue = cmhtest(arr)
        return np.log2(odds_ratio) if odds_ratio != 'NA' else 'NA', pvalue
    
    if n_sam == 2:
        # use fisher exact test
        idx_ppc1 = col_idx['ppc']['cond1'][0]
        idx_ppc2 = col_idx['ppc']['cond2'][0]
        idx_gbc1 = col_idx['gbc']['cond1'][0]
        idx_gbc2 = col_idx['gbc']['cond2'][0]
        data_out = data_pass.iloc[:, idx_gene_cols].copy()
        data_out['pvalue'] = data_pass.apply(lambda x: get_pvalue_2_sample(x, idx_ppc1, idx_ppc2, idx_gbc1, idx_gbc1), axis=1)
        data_out = add_FDR_col(data_out, 'pvalue')
        # logger.info(data_out.head())

        
        # odds_ratio = (ppc2/ppc1) / (gbc2/gbc1) = ppc2 * gbc1 / (ppc1 * gbc2)
        def get_odds_ratio(row):
            ppc1 = row[idx_ppc1]
            ppc2 = row[idx_ppc2]
            gbc1 = row[idx_gbc1]
            gbc2 = row[idx_gbc2]
            tmp = ppc1 * gbc2
            if all([ppc1, ppc2, gbc1, gbc2, tmp]):
                return ppc2 * gbc1 / (tmp)
            return np.nan

        try:
            # log2fc = data_pass.apply(lambda x: np.log2(x[idx_ppc2] * x[idx_gbc1] / (x[idx_ppc1] * x[idx_gbc2])), axis=1)
            log2fc = data_pass.apply(lambda x: np.log2(get_odds_ratio(x)), axis=1)
        except:
            e = traceback.format_exc()
            logger.warning(e)
            x = list(data_pass.iloc[0])
            logger.warning(x)
            logger.info(get_odds_ratio(x))
            logger.warning(f'{x[idx_ppc2]=}, {x[idx_gbc1]=}, {x[idx_ppc1]=}, {x[idx_gbc2]=}')
            sys.exit(1)
        
        data_out['log2fc'] = log2fc

    # run the cmhtest, the data should have 3 dimensions
    else:
        # logger.info(data_pass.columns)
        # logger.info(data_pass.head().values)
        data_out = data_pass.iloc[:, idx_gene_cols].copy()
        
        data_out[['log2fc', 'pvalue']] = data_pass.iloc[:, n_gene_cols:].apply(cmhtest_for_pindex_change, axis=1, result_type='expand')
        data_out = add_FDR_col(data_out, 'pvalue')
        data_out = data_out.sort_values('FDR')

    data_out_full = pd.concat([data_out, data_drop.iloc[:, idx_gene_cols]]).fillna('NA')
    data_out_full.to_csv(fno, sep='\t', index=False, na_rep='NA')
    # logger.info(f'pindex_change done : {ana_type}')
    # logger.info(data_out.head())

def add_value_to_gtf(gene_info, pro_up, pro_down, gb_down_distance, tts_padding, tts_down_length=5000, islongerna=False):
    strand = gene_info['strand']
    gene_raw_s, gene_raw_e = gene_info['start'], gene_info['end']
    if not islongerna:
        last_exon_s, last_exon_e = gene_info['last_exon']
    else:
        last_exon_s, last_exon_e = None, None
    if strand == '+':
        pp_start = gene_raw_s - pro_up
        pp_end = gene_raw_s + pro_down - 1
        gb_start = gene_raw_s + gb_down_distance
        gb_end = gene_raw_e
        tts = gene_raw_e
        strand_idx = 0
        tts_start = gene_raw_e
        tts_end = gene_raw_e + tts_padding
        tts_down_region_s, tts_down_region_e = gene_raw_e, gene_raw_e + tts_down_length
    else:
        pp_start = gene_raw_e - (pro_down - 1)
        pp_end = gene_raw_e + pro_up
        gb_start = gene_raw_s
        tts = gene_raw_s
        gb_end = gene_raw_e - gb_down_distance
        strand_idx = 1
        tts_start = gene_raw_s - tts_padding
        tts_end = gene_raw_s
        tts_down_region_s, tts_down_region_e = gene_raw_s - tts_down_length, gene_raw_s
    
    # modify here
    # skip the get gene_seq step
    gene_seq = None
    gb_seq_N = 0

    # gene_seq = get_seqence_from_fa(fa_idx, fh_fa, chr_, gene_raw_s, gene_raw_e)# already upper case
    # if gene_seq.count('N') == 0:
    #     gene_seq = None  # no need to check the mappable sites
    #     gb_seq_N = 0
    # else:
    #     gb_seq_N = gene_seq[gb_start - gene_raw_s: gb_end - gene_raw_s + 1].count('N')
    
    gb_len_mappable = gb_end - gb_start + 1 - gb_seq_N

    new_info = dict(zip(
        ['pp_start', 'pp_end', 'gb_start', 'gb_end', 'strand_idx', 'gb_len_mappable', 'gene_seq', 'tts_start', 'tts_end', 'last_exon_s', 'last_exon_e', 'tts_down_region_s', 'tts_down_region_e'], 
        [ pp_start,   pp_end,   gb_start ,  gb_end ,  strand_idx ,  gb_len_mappable,   gene_seq, tts_start, tts_end, last_exon_s, last_exon_e, tts_down_region_s, tts_down_region_e]
        ))
    gene_info.update(new_info)
    return gene_info

def build_tss(gtf_info, fn_tss, fn_tss_tts):
    # create the TSS file
    # chr1	11874	11874	NR_046018.2	+
    logger.info('creating TSS file')
    logger.debug(f'fn_tss = {fn_tss}, fn_tss_tts = {fn_tss_tts}')
    with open(fn_tss, 'w') as f, open(fn_tss_tts, 'w') as o2:
        for k, v in gtf_info.items():
            itss, itts = [v['start'], v['end']] if v['strand'] == '+' else [v['end'], v['start']]
            
            f.write(f'{v["chr"]}\t{itss}\t{itss}\t{k}\t{v["strand"]}\n')
            o2.write(f'{v["chr"]}\t{itss}\t{itts}\t{k}\t{v["strand"]}\n')

def process_gtf(home, fn_gtf, pwout, force_rebuild=False, fake_gtf_path=None):
    """
    process the gtf file, get the gene start and end position, and other info
    gtf position is 1-idx, full closed
    fake_gtf_path is just used for building the pre-built pkl file inside of the docker container
    """
    custom_temp = f"{pwout}/temp"
    logger.debug(f'set temp folder to {custom_temp}')
    os.makedirs(custom_temp, exist_ok=True)
    tempfile.tempdir = custom_temp  # This is a global setting
    os.environ['TMPDIR'] = custom_temp  # This is also global

    err = {'no_transcript_id': 0, 'no_gene_name': 0, 'invalid_line_format': 0, 'invalid_strand': 0}
    fn_gtf = os.path.realpath(fn_gtf)
    fn_gtf_lb = os.path.basename(fn_gtf).replace('.gz', '').replace('.gtf', '')
    
    if fake_gtf_path is not None:
        fn_gtf_md5 = get_md5_str(fake_gtf_path)[:8]
    else:
        fn_gtf_md5 = get_md5_str(fn_gtf)[:8]
    fn_gtf_size = os.path.getsize(fn_gtf)
    fn_gtf_lb = f'{fn_gtf_lb}_md5_{fn_gtf_md5}_size_{fn_gtf_size}'
    logger.debug(f'input gtf file = {fn_gtf}, file name string md5sum = {fn_gtf_md5}, file size = {fn_gtf_size}')
    
    logger.debug({'fn_gtf': fn_gtf, 'pwout': pwout, 'force_rebuild': force_rebuild, 'fn_gtf_lb': fn_gtf_lb})
    
    # logger.warning(fn_gtf)
    fn_tss = f'{pwout}/intermediate/{fn_gtf_lb}.tss.txt'
    fn_tss_tts = f'{pwout}/intermediate/{fn_gtf_lb}.tss_tts.txt'
    
    
    # check if have write permission to the folder of the gtf file
    if not os.access(os.path.dirname(fn_gtf), os.W_OK) or fn_gtf.startswith('/app/nrsa'):
        gtf_pwout = f'{home}/.nrsa'
        os.makedirs(gtf_pwout, exist_ok=True)
    else:
        gtf_pwout = os.path.dirname(fn_gtf)
    
    if 'gtf' not in fn_gtf.rsplit('.', 2)[-2:]:
        logger.error(f'the gtf file should have the extension of .gtf: {fn_gtf}')
        sys.exit(1)
    fn_gtf_pkl_default = f'{os.path.dirname(fn_gtf)}/{fn_gtf_lb}.gtf_info.pkl'
    fn_gtf_pkl = f'{gtf_pwout}/{fn_gtf_lb}.gtf_info.pkl'
    fn_gtf_meta_json = f'{gtf_pwout}/{fn_gtf_lb}.gtf_meta.json'
    
    if os.path.exists(fn_gtf_pkl_default):
        fn_gtf_pkl = fn_gtf_pkl_default
    logger.debug(f'fn_gtf_pkl = {fn_gtf_pkl}, fn_gtf_meta_json = {fn_gtf_meta_json}')
    if os.path.exists(fn_gtf_meta_json):
        # print the content
        with open(fn_gtf_meta_json) as f:
            logger.debug(f'gtf_meta_json = \n{f.read()}')

    if os.path.exists(fn_gtf_pkl) and not force_rebuild:
        logger.debug(f'loading gtf from pickle: {fn_gtf_pkl}')
        
        with open(fn_gtf_pkl, 'rb') as f:
            gtf_info = pickle.load(f)
        if pwout is not None and (not os.path.exists(fn_tss) or not os.path.exists(fn_tss_tts)):
            build_tss(gtf_info, fn_tss, fn_tss_tts)

        return gtf_info, fn_tss, fn_tss_tts, err
    
    gtf_col_idx = {
        'chr': 0,
        'source': 1,
        'feature': 2,
        'start': 3,
        'end': 4,
        'score': 5,
        'strand': 6,
        'attribute': 8,
    }
    res_raw = {}  # k1 = gene name, k2 = transcript_id, v = {'chr': '', 'strand': '', 'gene_id': '', 'gene_name': '', 'start': 0, 'end': 0}
    
    last_exon = {}
    with gzip.open(fn_gtf, 'rt') if fn_gtf.endswith('.gz') else open(fn_gtf) as f:
        while True:
            i = f.readline()
            if i[0] != '#' or not i:
                break
        f.seek(f.tell() - len(i))
        skipped_transcript_on_random_chr = 0
        for i in f:
            line = i.strip().split('\t')
            # chr1    hg19_ncbiRefSeq    exon    66999252    66999355    0.000000    +    .    gene_id "NM_001308203.1"; transcript_id "NM_001308203.1"; gene_name "SGIP1";
            # 1       ensembl_havana  gene    11869   14412   .       +       .       gene_id "ENSG00000223972"; gene_version "4"; gene_name "DDX11L1"; gene_source "ensembl_havana"; gene_biotype "pseudogene";
            line_err = None
            try:
                chr_, region_type, start, end, strand, attribute_raw = [line[_] for _ in [gtf_col_idx[_] for _ in ['chr', 'feature', 'start', 'end', 'strand', 'attribute']]]
                start = int(start)
                end = int(end)
            except:
                line_err = 'invalid_line_format'
            if region_type != 'exon':
                continue
            if strand not in {'+', '-'}:
                line_err = 'invalid_strand'
            chr_ = refine_chr(chr_)

            m = re.search(r'transcript_id\s+"?([\w._\-]+)', attribute_raw)
            if m:
                transcript_id = m.group(1)
            else:
                line_err = 'no_transcript_id'
                
            m = re.search(r'gene_name\s+"?([\w._\-]+)', attribute_raw)
            if m:
                gene_name = m.group(1)
            else:
                line_err = 'gene_name_not_found'
            
            if line_err:
                tmp = err.setdefault(line_err, 0)
                err[line_err] += 1
                if err[line_err] < 10:
                    logger.debug(f'error during processing gtf {line_err}: {i}')
                
                continue
            
            if '_' in chr_:
                skipped_transcript_on_random_chr += 1
                continue # skip the line if the chr_ contains underscore (e.g. chr1_KI270706v1_random)

            # {'chr': '', 'strand': '', 'gene_name': '', 'start': 0, 'end': 0}
            ires_exon = {'chr': chr_, 'strand': strand, 'gene_name': gene_name, 'start': start, 'end': end}
            ires = res_raw.setdefault(gene_name, {}).setdefault(transcript_id, ires_exon)
            if strand == '+':
                # if plus strand, will keep changing, leaving only the last appearance
                last_exon[transcript_id] = (start, end)
            elif transcript_id not in last_exon:
                # if on minus strand, will only try if transcript id is not found before
                last_exon[transcript_id] = (start, end)
                
            
            # deal with the case when the same transcript-ID with different chr or strand
            if chr_ != ires['chr'] or strand != ires['strand']:
                transcript_id_new = f'{transcript_id}@{chr_}@{strand}'
                ires = res_raw.setdefault(gene_name, {}).setdefault(transcript_id_new, ires_exon)
                if strand == '+':
                    # if plus strand, will keep changing, leaving only the last appearance
                    last_exon[transcript_id_new] = (start, end)
                elif transcript_id_new not in last_exon:
                    # if on minus strand, will only try if transcript id is not found before
                    last_exon[transcript_id_new] = (start, end)

            if start < ires['start']:
                ires['start'] = start
            if end > ires['end']:
                ires['end'] = end

    # merge transcripts with same start and end position
    res = {}
    n_merged = 0
    meta = {'initial_n_transcripts': 0, 'n_merged': 0, 'final_n_transcripts': 0, 'skipped_transcript_on_random_chr': skipped_transcript_on_random_chr}

    for gn, v1 in res_raw.items():
        tmp = {}  # key = unique_id
        for transcript_id, v2 in v1.items():
            start, end, strand = v2['start'], v2['end'], v2['strand']
            unique_id = f'{v2["chr"]}_{start}_{end}_{strand}'
            v2['tss'] = start if strand == '+' else end
            v2['tts'] = end if strand == '+' else start
            tmp.setdefault(unique_id, []).append(transcript_id) # these transcripts have the same start and end position
            meta['initial_n_transcripts'] += 1
        for transcript_list in tmp.values():
            if len(transcript_list) > 1:
                n_merged += len(transcript_list) - 1
                
            transcript_id_new = ';'.join(transcript_list)
            transcript_id_demo = transcript_list[0]
            res[transcript_id_new] = v1[transcript_id_demo]
            res[transcript_id_new]['last_exon'] = last_exon[transcript_id_demo]
    meta['n_merged'] = n_merged
    meta['final_n_transcripts'] = len(res)
    logger.debug(f'g@merged {n_merged} transcripts with same start and end position')
    # logger.debug(meta)
    
    with open(fn_gtf_pkl, 'wb') as o:
        logger.debug(f'dump to {fn_gtf_pkl}')
        pickle.dump(res, o)
    
    with open(fn_gtf_meta_json, 'w') as o:
        logger.debug(f'dump to {fn_gtf_meta_json}')
        json.dump(meta, o, indent=4)

    if pwout is not None:
        build_tss(res, fn_tss, fn_tss_tts)
    err_total = sum(err.values())
    if err_total:
        logger.info(f'error in parsing gtf file: {err}')
    
    return res, fn_tss, fn_tss_tts, err



def check_dependency():
    """
    check if the dependency: HOMER and bedtools are installed
    """
    # check HOMER
    err = 0

    retcode = run_shell('which makeTagDirectory')
    if retcode != 0:
        logger.error("HOMER is not installed, please install HOMER first")
        err = 1

    # check bedtools
    retcode = run_shell('which bedtools')
    if retcode != 0:
        logger.error("bedtools is not installed, please install bedtools first")
        err = 1
    
    # check python packages
    required_packages = ['pydeseq2', 'pandas', 'numpy', 'statsmodels', 'scipy', 'matplotlib', 'fisher']
    import importlib.util
    missing_python_packages = []
    for pkg in required_packages:
        if importlib.util.find_spec(pkg) is None:
            missing_python_packages.append(pkg)
            err = 1
    if missing_python_packages:
        logger.error(f"Missing python packages: {missing_python_packages}, please install them first.")
    
    return err


def run_shell(cmd, echo=False, ret_info=False):
    """
    run shell command and check the return code and stdout stderr
    """
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    s = time.time()
    stdout, stderr = p.communicate()
    duration = time.time() - s
    retcode = p.returncode
    if duration > 2:
        logger.debug(f'cmd = {cmd}, retcode = {retcode}, duration = {duration:.2f}s')
    if ret_info:
        echo = False
    if retcode:
        logger.error(f"Error running command: {cmd}")
        logger.error(f"stdout: {stdout.decode()}")
        logger.error(f"stderr: {stderr.decode()}")
    if not ret_info:
        return retcode
    return retcode, stdout, stderr

def get_lb(fn):
    lb = os.path.basename(fn).replace('.gz', '').rsplit('.', 1)[0]
    return re.sub(r'[\._\-]*sort(ed)?\b', '', lb)



def pre_count_for_bed_legacy(fn_lb, fn_out_bed, pw_bed, bin_size=200, reuse=True):
    """
    process the bed file, and get the count of the read end regarding strand.
    will export 2 files, 
    1. dict, key = chr, k2 = strand, k3 = chunk_number. value = count of reads with read end in this chunk.  the chunk_number is got from position // bin_size  (default chunksize = 10bp)
    2. bed file, record the read end chr, read_end position, strand and the count of reads with this read end position
    
    bin_size = 200 # compared bin size of 10, 20, 50, 100 and 200 and 500.  200 is the best, maybe because when getting the gb count, the size fit the distribution
    return chr_map is used to make the chr pattern in the gtf file and input bed file consistant, esp. when running bedtools

    """
    fn_count_bin = f'{pw_bed}/{fn_lb}.count.bin_of_{bin_size}.pkl'
    fn_count_per_base = f'{pw_bed}/{fn_lb}.count.per_base.pkl'
    fn_chr_map = f'{pw_bed}/{fn_lb}.chr_map.pkl'
    fn_n_lines = f'{pw_bed}/{fn_lb}.line_count.txt'
    
    if reuse and os.path.exists(fn_count_per_base) and os.path.getsize(fn_count_per_base) > 10 and os.path.exists(fn_count_bin) and os.path.getsize(fn_count_bin) > 10:
        logger.info(f'Loading pre-counting data - {fn_lb}...')
        with open(fn_count_per_base, 'rb') as f:
            count_per_base = pickle.load(f)
        with open(fn_count_bin, 'rb') as f:
            count_bin = pickle.load(f)
        logger.info('Loaded.')
        return count_per_base, count_bin

    count_bin = {'bin_size': bin_size}
    count_per_base = {}
    logger.info(f'Pre-counting for {fn_lb}')
    n_lines = 0
    n_error = 0
    single_col = 0
    fatal = 0
    with gzip.open(fn_out_bed, 'rt') if fn_out_bed.endswith('.gz') else open(fn_out_bed, 'r') as fh_bed:
        for i in fh_bed:
            # chr1	10511	10569	A00758:60:HNYK7DSXX:4:2461:27181:23171	32	-
            try:
                chr_, s, e, _, _, strand, *_ = i[:-1].split('\t', 6)
            except:
                nchr, ncol = len(i), len(i.split('\t'))
                logger.debug(f'Invalid line format: nchar = {nchr}, ncol = {ncol}, first 100 char = {i[:100]}')
                if ncol == 1:
                    single_col += 1
                else:
                    n_error += 1
                    if n_error > 100:
                        fatal = 1
                        break
                continue
            # because we are examining where the transcript stops/pauses
            # so we only care about the read end, no matter if we will use it in TSS or TTS or gene body calculation
            strand_idx, read_end = (0, int(e)) if strand == '+' else (1, int(s) + 1)

            chunk = read_end // bin_size
            count_bin.setdefault(chr_, {}).setdefault(chunk, [0, 0])
            count_bin[chr_][chunk][strand_idx] += 1

            count_per_base.setdefault(chr_, {}).setdefault(read_end, [0, 0])  # for each read end, the count of + and - strand
            count_per_base[chr_][read_end][strand_idx] += 1
            n_lines += 1

        # refine the chr_
        chr_map = {refine_chr(k): k for k in count_bin}
        count_per_base = {refine_chr(k): v for k, v in count_per_base.items()}
        count_bin = {refine_chr(k): v for k, v in count_bin.items()}
    if n_error:
        prefix = 'more than ' if fatal else ''
        logger.error(f'{prefix} {n_error} lines with invalid format were found, check the input file: {fn_out_bed} and \nlog file for details.\nthe first 6 columns should be chr, start, end, read_name, score, strand')
        logger.error(f'now exiting...')
        sys.exit(1)
    logger.info('done.')
    # dump the pickle
    
    
    logger.debug(f'dumping to pickle: {fn_count_per_base}, {fn_count_bin}, {fn_chr_map}')
    with open(fn_count_per_base, 'wb') as f:
        pickle.dump(count_per_base, f)
    with open(fn_count_bin, 'wb') as f:
        pickle.dump(count_bin, f)
        
    with open(fn_chr_map, 'wb') as f:
        pickle.dump(chr_map, f)
    logger.debug('dump to pickle done.')
    
    with open(fn_n_lines, 'w') as f:
        f.write(f'{n_lines}\n')
    
    
    return count_per_base, count_bin


def pre_count_for_bed(fn_lb, fn_out_bed, pw_bed, bin_size=200, reuse=True):
    """
    process the bed file, and get the count of the read end regarding strand.
    will export 2 files, 
    1. dict, key = chr, k2 = strand, k3 = chunk_number. value = count of reads with read end in this chunk.  the chunk_number is got from position // bin_size  (default chunksize = 10bp)
    2. bed file, record the read end chr, read_end position, strand and the count of reads with this read end position
    
    bin_size = 200 # compared bin size of 10, 20, 50, 100 and 200 and 500.  200 is the best, maybe because when getting the gb count, the size fit the distribution
    return chr_map is used to make the chr pattern in the gtf file and input bed file consistant, esp. when running bedtools

    """
    fn_chr_map = f'{pw_bed}/{fn_lb}.chr_map.pkl'
    fn_n_lines = f'{pw_bed}/{fn_lb}.line_count.txt'
    
    fn_pre_count_list_json = f'{pw_bed}/{fn_lb}.pre_count_list.json'
    
    if reuse and os.path.exists(fn_pre_count_list_json) and os.path.getsize(fn_pre_count_list_json) > 10:
        ok = 1
        with open(fn_pre_count_list_json) as f:
            file_list = json.load(f)
            for chr_, (fn_count_per_base, fn_count_bin) in file_list.items():
                if not os.path.exists(fn_count_per_base) or not os.path.exists(fn_count_bin):
                    ok = 0
                    break
        if ok:
            logger.info(f'Reuse existing pre-counting results - {fn_lb}...')
            return file_list

    file_list = {} # key = chr, v = fn_count_per_base, fn_count_bin
    chr_map = {}
    logger.info(f'Pre-counting for {fn_lb}')
    n_lines = 0
    n_error = 0
    single_col = 0
    fatal = 0
    count_bin = {'bin_size': bin_size}
    count_per_base = {}
    chr_prev = None
    
    
    def dump_res(chr_prev, count_per_base, count_bin):
        # dump the pickle
        chr_prev_new = refine_chr(chr_prev)
        chr_map[chr_prev_new] = chr_prev
        fn_count_bin = f'{pw_bed}/{fn_lb}.count.bin_of_{bin_size}_chr{chr_prev_new}.pkl'
        fn_count_per_base = f'{pw_bed}/{fn_lb}.count.per_base_chr{chr_prev_new}.pkl'
        logger.debug(f'dumping to pickle: {fn_count_per_base}, {fn_count_bin}')
        with open(fn_count_per_base, 'wb') as f:
            pickle.dump(count_per_base, f)
        with open(fn_count_bin, 'wb') as f:
            pickle.dump(count_bin, f)
        file_list[chr_prev_new] = (fn_count_per_base, fn_count_bin)
        logger.debug('dump to pickle done.')
        count_bin = {'bin_size': bin_size}
        count_per_base = {}
        return count_per_base, count_bin
    
    
    with gzip.open(fn_out_bed, 'rt') if fn_out_bed.endswith('.gz') else open(fn_out_bed, 'r') as fh_bed:
        for i in fh_bed:
            # chr1	10511	10569	A00758:60:HNYK7DSXX:4:2461:27181:23171	32	-
            try:
                chr_, s, e, _, _, strand, *_ = i[:-1].split('\t', 6)
            except:
                nchr, ncol = len(i), len(i.split('\t'))
                logger.debug(f'Invalid line format: nchar = {nchr}, ncol = {ncol}, first 100 char = {i[:100]}')
                if ncol == 1:
                    single_col += 1
                else:
                    n_error += 1
                    if n_error > 100:
                        fatal = 1
                        break
                continue
            # because we are examining where the transcript stops/pauses
            # so we only care about the read end, no matter if we will use it in TSS or TTS or gene body calculation
            strand_idx, read_end = (0, int(e)) if strand == '+' else (1, int(s) + 1)
            if chr_prev and chr_ != chr_prev:
                count_per_base, count_bin = dump_res(chr_prev, count_per_base, count_bin)


            chunk = read_end // bin_size
            count_bin.setdefault(chunk, [0, 0])
            count_bin[chunk][strand_idx] += 1

            count_per_base.setdefault(read_end, [0, 0])  # for each read end, the count of + and - strand
            count_per_base[read_end][strand_idx] += 1
            n_lines += 1
            chr_prev = chr_
    # dump the last chr
    if chr_prev:
        count_per_base, count_bin = dump_res(chr_prev, count_per_base, count_bin)

    if n_error:
        prefix = 'more than ' if fatal else ''
        logger.error(f'{prefix} {n_error} lines with invalid format were found, check the input file: {fn_out_bed} and \nlog file for details.\nthe first 6 columns should be chr, start, end, read_name, score, strand')
        logger.error(f'now exiting...')
        sys.exit(1)
    logger.info('done.')

    with open(fn_chr_map, 'wb') as f:
        pickle.dump(chr_map, f)
    logger.debug('dump to pickle done.')
    
    with open(fn_n_lines, 'w') as f:
        f.write(f'{n_lines}\n')
    with open(fn_pre_count_list_json, 'w') as o:
        json.dump(file_list, o, indent=3)
    return file_list


def process_input(pwout_raw, fls, respect_sorted=False):
    """
    process input files, convert bam to bed if needed, and build the index for bed files
    return a list of [fn_bed, idx_bed, fn_lb]
    """

    if fls is None or not fls:
        return None
    fls_new = []
    err = 0
    for fn in fls:
        fn = os.path.realpath(fn)
        if not os.path.exists(fn):
            logger.error(f'Input file not exist: {fn}')
            err = 1
            continue
        fls_new.append(fn)
    if err:
        sys.exit(1)
    fls = fls_new
    if err:
        sys.exit(1)
    fls = fls_new
    
    err = 0
    res = []
    lb_map = {}
    if not os.path.exists(f'{pwout_raw}/bed'):
        os.makedirs(f'{pwout_raw}/bed', exist_ok=True)
    for fn in fls:
        fn_lb = get_lb(fn)
        gz_suffix = '.gz' if fn.endswith('.gz') else ''
        fn_out_bed = f'{pwout_raw}/bed/{fn_lb}.sorted.bed'
        fn_out_bed_gz = f'{pwout_raw}/bed/{fn_lb}.sorted.bed{gz_suffix}'

        if os.path.exists(fn_out_bed) and os.path.getsize(fn_out_bed) > 10:
            logger.debug(f'bed file already exist: {fn_out_bed}')
            res.append([fn_lb, fn_out_bed])
            continue
        elif gz_suffix and os.path.exists(fn_out_bed_gz) and os.path.getsize(fn_out_bed_gz) > 10:
            res.append([fn_lb, fn_out_bed_gz])
            logger.debug(f'bed file already exist (gzip): {fn_out_bed_gz}')
            continue
        logger.debug(f'bed file not exist yet: {fn_out_bed}')
        
        fn_for_check = re.sub(r'\.gz$', '', fn)
        if fn_for_check.endswith('.bam'):
            file_size = os.path.getsize(fn) / 1024 / 1024/1024 # GB
            logger.info(f'Converting bam to bed: {fn}, bam file size = {file_size:.1f} GB')
            # still need to do the sorting, because for the following steps using bedtools coverage, the sorted bed will be more memory efficient
            # bedtools sort is faster than linux sort
            fn_converted_bed = f'{pwout_raw}/bed/{fn_lb}.converted.bed'
            cmd = f"""bedtools bamtobed -i {fn} > {fn_converted_bed}"""
            retcode = run_shell(cmd)
            if not os.path.exists(fn_out_bed) or retcode:
                logger.error(f'Error converting bam to bed: {fn}')
                err = 1
            else:
                file_size = os.path.getsize(fn_converted_bed) / 1024 / 1024/1024 # GB
                logger.debug(f'converted bed file size: {file_size:.2f} GB')
            fn_for_check = fn_converted_bed

        # wrap the converted bed too
        if fn_for_check.endswith('.bed'):
            # check if sorted
            if respect_sorted and 'sort' in fn:
                is_sorted = 1
            else:
                is_sorted = check_is_sorted(fn)
            if is_sorted:
                fn_abs = os.path.realpath(fn)
                fn_dest = fn_out_bed_gz if gz_suffix else fn_out_bed
                logger.debug(f'Creating symlink for sorted bed file: {fn_abs} -> {fn_dest}')
                force_symlink(fn_abs, fn_dest)
                ires = [fn_lb, fn_dest]
            else:
                logger.info(f'Sorting input file: {fn}')
                # bedtools can handle gzip format
                cmd = f'bedtools sort -i {fn} > {fn_out_bed}'
                retcode = run_shell(cmd)
                if retcode:
                    logger.error(f'Error sorting bed file: {fn}')
                    err = 1
                    continue
                ires = [fn_lb, fn_out_bed]
        else:
            logger.error(f"Input file '{fn}' should be in bed or bam format")
            err = 1
            continue
        res.append(ires)
    if err:
        return None
    return res



def get_overlap(s1, e1, s2, e2):
    """
    get the overlap fraction of two intervals
    """
    if s1 > e2 or s2 > e1:
        return 0
    return min(e1, e2) - max(s1, s2) + 1

def gtf_compare(gtf_info, fn_peak_gtf):
    """
    verified, results match with perl code
    compare the refseq_gtf and peak_gtf, and get the overlap fraction
    overlap_frac, removed.
    """
    # gtf_info, will use the start and end

    # fn_peak_gtf is like below
    # chr1	homer	exon	4768429	4785693	0.622000	-	.	gene_id "chr1-1854-1"; transcript_id "chr1-1854-1"
    # chr1	homer	exon	4785562	4788480	0.040000	+	.	gene_id "chr1-2-0"; transcript_id "chr1-2-0"
    ts_new_all = {}
    ts_by_chr = {}
    with open(fn_peak_gtf) as f:
        for i in f:
            line = i[:-1].split('\t')
            chr_, _, _, start, end, _, strand = line[:7]
            chr_ = refine_chr(chr_)
            # line[-1] is like gene_id "chr1-3785-1"; transcript_id "chr1-3785-1"
            m = re.search(r'transcript_id "(\S+)"', line[-1])
            if m:
                transcritp_id = m.group(1)
            else:
                continue
            chr_strand_key = f'{chr_}@{strand}'
            ts_new_all[transcritp_id] = {'chr': chr_, 'start': int(start), 'end': int(end), 'strand': strand}
            ts_by_chr.setdefault(chr_strand_key, set()).add(transcritp_id)

    # annotated genes
    # ref_gene_type = {'coding': [], 'noncoding': []}
    ts_new_with_overlap = set()
    for gene_info in gtf_info.values():
        ref_chr, ref_strand, ref_s, ref_e = gene_info['chr'], gene_info['strand'], gene_info['start'], gene_info['end']
        # overlap_sum = 0
        # ref_gene_len = ref_e - ref_s + 1
        chr_strand_key = f'{ref_chr}@{ref_strand}'
        ts_list = ts_by_chr.get(chr_strand_key, set())
        for its in ts_list:
            its_info = ts_new_all[its]
            overlap = get_overlap(ref_s, ref_e, its_info['start'], its_info['end'])
            # overlap_sum += overlap
            if overlap > 0:
                ts_new_with_overlap.add(its)
            
        # if overlap_sum >= overlap_frac * ref_gene_len:
        #     itype = 'noncoding' if ts.startswith('NR') else 'coding'
        #     ref_gene_type[itype].append(ts)

    # antisense genes, check overlap by flip the strand of the refseq transcripts and check overlap
    ts_new_no_overlap = set(ts_new_all) - ts_new_with_overlap # only cehck antisense for these transcripts
    # ts_new_antisense_overlap = set()
    for gene_info in gtf_info.values():
        ref_chr, ref_strand, ref_s, ref_e = gene_info['chr'], gene_info['strand'], gene_info['start'], gene_info['end']
        # ref_gene_len = ref_e - ref_s + 1
        ref_strand_flip = '-' if ref_strand == '+' else '+'
        chr_strand_key = f'{ref_chr}@{ref_strand_flip}'
        ts_list = ts_by_chr.get(chr_strand_key, set()) - ts_new_with_overlap
        # overlap_sum = 0
        # its_l = set()
        for its in ts_list:
            its_info = ts_new_all[its]
            overlap = get_overlap(ref_s, ref_e, its_info['start'], its_info['end'])
            # overlap_sum += overlap
            # its_l.add(its)
            if overlap > 0:
                ts_new_with_overlap.add(its)
    
        # if overlap_sum >= overlap_frac * ref_gene_len:
        #     ts_new_antisense_overlap |= its_l

    # divergent genes
    ts_new_no_overlap -= ts_new_with_overlap
    for gene_info in gtf_info.values():
        ref_chr, ref_strand, ref_s, ref_e = gene_info['chr'], gene_info['strand'], gene_info['start'], gene_info['end']
        ref_gene_len = ref_e - ref_s + 1
        ref_strand_flip = '-' if ref_strand == '+' else '+'
        if ref_gene_len < 1000:
            continue
        if ref_strand == '+':
            ref_strand_flip = '-'
            ref_s_new = ref_s -500
            ref_e_new = ref_s + 499
        else:
            ref_strand_flip = '+'
            ref_s_new = ref_e - 499
            ref_e_new = ref_e + 500
        # overlap_sum = 0
        # trans_len = 0
        # its_l = set()
        chr_strand_key = f'{ref_chr}@{ref_strand_flip}'
        ts_list = ts_by_chr.get(chr_strand_key, set()) - ts_new_with_overlap
        for its in ts_list:
            its_info = ts_new_all[its]
            overlap = get_overlap(ref_s_new, ref_e_new, its_info['start'], its_info['end'])
            # overlap_sum += overlap
            # trans_len += its_info['end'] - its_info['start'] + 1
            # its_l.add(its)
            if overlap > 0:
                ts_new_with_overlap.add(its)
        # if overlap_sum >= 0.1 * trans_len and overlap_sum <= 0.5 * ref_gene_len:
        #     divergent.setdefault(ts, []).append(its)
    ts_new_no_overlap -= ts_new_with_overlap
    
    # with open('other_genes.v2.txt', 'w') as o:
    #     print('\n'.join(sorted(ts_new_no_overlap)), file=o)
    
    return ts_new_no_overlap

def get_other_region(other_genes, fn_peak_txt):
    """
    verified, results are same as perl output
    """
    # other_genes is ts_new_no_overlap above , is a set
    # #PeakID chr     start   end     strand  Initial read depth      length (bp)
    # chr16-266-1     chr16   57390979        57391720        -       58.1    741.000
    # chr5-1675-0     chr5    146261022       146262398       +       49.2    1376.000
    
    regions = []
    # the chr_ here is refined format
    with open(fn_peak_txt) as f:
        for i in f:
            if i[0] == '#':
                continue
            line = i[:-1].split('\t')
            if line[0] in other_genes:
                line[1] = refine_chr(line[1])  # line[1] is chr
                line[2], line[3] = int(line[2]), int(line[3])
                regions.append(line[1:5])
    logger.debug(f'input other_genes = {len(other_genes)}, regions added info = {len(regions)}')
    return regions


def is_in_any_interval(intervals, v):
    """
    check if a specific number fall in any interval
    the intervals must be sorted
    inf_neg is the negative infinity, defined outside
    """
    
    i = bisect.bisect_left(intervals, [v, inf_neg])
    if i < len(intervals) and intervals[i][0] <= v <= intervals[i][1]:
        return True
    if i > 0 and intervals[-1][0] <= v <= intervals[-1][1]:
        return True
    return False


def get_enhancer(other_region, fn_fantom, fn_association, fn_tss_tts, lcut=400, thres_long_eRNA=10_000, distance_for_merge=500, filter=False):
    """
    equiv to the central function in eRNA.pl
    the chr_ here in other_regions and the files are original format, don't need to be refiend
    thres_long_eRNA, only the length > this value, will be considered as longeRNA
    lcut, cutoff of 5' distance
    distance_for_merge: if 2 eRNA regions distance < 500, merge these 2
    filter:  filter the result based on TSS-TTS
    chr_ here is refined
    """
    # other_region is a list of [chr, start, end, strand], got from above get_other_region function
    regions = {}
    for i in other_region:
        chr_, s, e, strand = i[:4]
        regions.setdefault(chr_, {'+': [], '-': []})[strand].append([s, e])

    fantom_sites= {} # k1 = chr, k2 = each pos
    enh_sites = {}

    logger.debug(f'{fn_fantom}, {fn_association}')
    if fn_fantom:
        # chr1	839741	840250	chr1:839741-840250	24	.	839787	839788	0,0,0	2	20,436	0,73
        # chr1	840753	841210	chr1:840753-841210	32	.	840811	840812	0,0,0	2	8,347	0,110
        with open(fn_fantom) as f:
            for i in f:
                chr_, s, e = i.split('\t')[:3]
                s, e = int(s), int(e)
                chr_ = refine_chr(chr_)
                if chr_ not in fantom_sites:
                    fantom_sites[chr_] = set()
                fantom_sites[chr_] |= set(range(s, e + 1))

    if fn_association:
        # actually the strand column in all rows are dot
        # 1:	#chrom	chr1
        # 2:	chromStart	66797292
        # 3:	chromEnd	67198741
        # 4:	name	chr1:67198280-67198800;NM_001037339;PDE4B;R:0.385;FDR:0
        # 5:	score	385

        with open(fn_association) as f:
            f.readline()
            for i in f:
                chr_, _, _, info = i[:-1].split('\t', 4)[:4]
                info = info.split(';')
                if len(info) < 5:
                    continue
                e_chr, e_start, e_end = re.split(r'[-:]', info[0])
                e_chr, e_start, e_end = refine_chr(e_chr), int(e_start), int(e_end)

                assoc_gn = info[2]
                if e_chr not in enh_sites:
                    enh_sites[e_chr] = {}
                for i in range(e_start, e_end + 1):
                    enh_sites[e_chr].setdefault(i, set()).add(assoc_gn)

        logger.debug(f'enh_sites chr keys = {sorted(enh_sites)}')
    
        # the gene list value to string 
        tmp = {}
        for k1, v1 in enh_sites.items():
            tmp[k1] = {}
            for k2, v2 in v1.items():
                tmp[k1][k2] = ','.join(sorted(v2))
        enh_sites = tmp
    

    # find center of enhancers
    enhancer_region = {}
    lerna_out = set()

    for chr_, v1 in regions.items():
        plus_region = sorted(v1['+'])
        minus_region = sorted(v1['-'], key=lambda x: x[1])
        minus_region_e = [_[1] for _ in minus_region] # used for jump to most adjacent plus region
        
        enhancer_region[chr_] = {} # k = e_plus, v = the e_minus

        for s_plus, e_plus in plus_region:
            len_plus_region = e_plus - s_plus
            is_long_eRNA_plus = True if len_plus_region > thres_long_eRNA else False
            iregion = {'start': e_plus, 'end': e_plus, 'center_list': []}
            
            # get the loop start for minus_region
            idx_minus_region = bisect.bisect_left(minus_region_e, s_plus - lcut) # from this index and after, all with d < lcut
            minus_region_sort_by_start = sorted(minus_region[idx_minus_region:])
            for s_minus, e_minus in minus_region_sort_by_start:
                if s_plus < s_minus and e_plus <= e_minus:
                    continue
                if s_minus > e_plus:
                    # -------->
                    #             <----------
                    break  
                
                distance = s_plus - e_minus
                if distance > lcut:
                    continue
                
                if s_minus < s_plus and e_minus <= e_plus:
                    # minus is before plus or have partial overlap
                    # ------------------->
                    #        <--------------------
                    center_pos = s_plus - (int(distance / 2))
                elif s_minus >= s_plus and e_minus < e_plus:
                    # minus region is inside of plus
                    center_pos = e_minus
                elif s_minus < s_plus and e_minus > e_plus:
                    # plus region is inside of minus
                    center_pos = s_plus
                else:
                    continue
                if s_minus < iregion['start']:
                    iregion['start'] = s_minus
                iregion['center_list'].append(str(center_pos))
                len_minus_region = e_minus - s_minus
                if len_minus_region > thres_long_eRNA:
                    lerna_out.add('\t'.join(map(str, [chr_, s_minus, e_minus, '-'])))
            # plus region is longeRNA
            n_centers = len(iregion['center_list'])
            if is_long_eRNA_plus and n_centers > 0:
                lerna_out.add('\t'.join(map(str, [chr_, s_plus, e_plus, '+'])))
            if n_centers > 0:
                # thre are matched minus strand to make a region
                if e_plus not in enhancer_region[chr_]:
                    iregion['center_list'] = sorted(iregion['center_list'])
                    enhancer_region[chr_][e_plus] = iregion
                else:
                    # merge them
                    if iregion['start'] < enhancer_region[chr_][e_plus]['start']:
                        enhancer_region[chr_][e_plus]['start'] = iregion['start']
                    enhancer_region[chr_][e_plus]['center_list'] = sorted(set(enhancer_region[chr_][e_plus]['center_list']) | set(iregion['center_list']))
    
    
    # combine adjacent regions
    enh_out = []
    n_merged_region = 0
    for chr_, v1 in enhancer_region.items():
        enh_out_chr = []
        for region_end in sorted(v1):
            iregion = v1[region_end]
            iregion['chr'] = chr_
            prev_region = enh_out_chr[-1] if enh_out_chr else None # last region
            if prev_region and iregion['start'] - prev_region['end'] < distance_for_merge:
                prev_region['end'] = iregion['end']
                prev_region['center_list'] = sorted(iregion['center_list'] + prev_region['center_list'])
                n_merged_region += 1
            else:
                enh_out_chr.append(iregion)
        
        for iregion in enh_out_chr:
            iregion['fantom_list'] = []
            iregion['asso_list'] = []
            for center_pos in iregion['center_list']:
                fantom_v, asso_v = 'N', 'NA'
                center_pos = int(center_pos)
                fantom_v = 'Y' if fantom_sites and chr_ in fantom_sites and center_pos in fantom_sites[chr_] else 'N'
                asso_v = 'NA'
                if enh_sites and chr_ in enh_sites and center_pos in enh_sites[chr_]:
                    asso_v = enh_sites[chr_][center_pos]
                
                iregion['fantom_list'].append(fantom_v)
                iregion['asso_list'].append(asso_v)
            enh_out.append(iregion)

    enh_out_str = ['\t'.join([i['chr'], str(i['start']), str(i['end']), ','.join(i['center_list']), ','.join(i['fantom_list']), ';'.join(i['asso_list'])]) for i in enh_out]
    lerna_out = sorted([_.split('\t') for _ in lerna_out], key=lambda x: (x[0], int(x[1]), int(x[2])))

    
    if filter:
        if not os.path.exists(fn_tss_tts):
            logger.warning(f'no tss_tts file found {fn_tss_tts}, skip filtering')
        else:
            tss_tts_info = process_tss_tts(fn_tss_tts)
            enh_out_str = [line_str for line_str in enh_out_str if filter_by_tss_tts(line_str.split('\t', 3), tss_tts_info)]
            lerna_out = [line for line in lerna_out if filter_by_tss_tts(line, tss_tts_info)]

    n_enhancer_raw = sum([len(v) for v in enhancer_region.values()])
    n_enhancer_after_merge = len(enh_out)
    exp_after_merge = n_enhancer_raw-n_merged_region
    n_enhancer_after_filtering = len(enh_out_str)
    logger.debug(f'enhancer_region, n_enhancer_raw = {n_enhancer_raw}, merged = {n_merged_region}, exp after_merge = {exp_after_merge},actual ={n_enhancer_after_merge}, match = {exp_after_merge == n_enhancer_after_merge}, n_enhancer_after_filtering={n_enhancer_after_filtering}')


    # enh_out = [chr, start, end, center_list, fantom_list, asso_list]
    return enh_out_str, lerna_out, 
            

def process_tss_tts(fn):
    """
    example, /Users/files/work/jb/work/NRSA_v2/new/ref/hg19/RefSeq-hg19-tss-tts.txt
    chr1	11874	14409	NR_046018.2	+
    chr1	3672278	3199733	XM_006495550.3	-
    this col2 can be bigger than col3, depending on the strand, not a typical bed format
    """
    tss_tts_info = {}
    with open(fn, 'r') as file:
        for line in file:
            temp = line.strip().split('\t')
            chr_ = refine_chr(temp[0])
            tss_tts_info.setdefault(chr_, {})[temp[3]] = {
                'chr': temp[0],
                'tss': int(temp[1]),
                'tts': int(temp[2]),
                'strand': temp[4]
            }

    return tss_tts_info


def filter_by_tss_tts(line, tss_tts_info, filter_tss=2000, filter_tts=20000):
    """
    if the enhancer is near any gene defined in tss_tts_info, then drop it
    equiv to the filter function in eRNA.pl
    return 1 if keep this line, 0 if drop
    line is a list like [chr, start, end, center_str, fantom_str, assocation_str]
    center_str and fantom_str are multiple elements joined by ',', assocation_str is multiple elements joined by ';'
    tss_tts_info is got by parsing tss_tts.txt file, a dict, contains chr, tss, tts, strand
    filter_tss: the minimum distance from enhancer to TSS
    filter_tts: the minimum distance from enhancer to TTS
    for the tss_tts_info, the tts can be smaller than tss, depending on the strand
    """
    # chr1	3361552	3377812	XR_865166.2	+
    # chr1	3672278	3199733	XM_006495550.3	-
    chr_, start, end = line[:3]
    start, end = int(start), int(end)
    if chr_ not in tss_tts_info:
        return 1
    keep = 1
    for v in tss_tts_info[chr_].values():
        if v['strand'] == '+':
            if start <= v['tts'] and (v['tss'] - end) < filter_tss:
                # end pos is within 2k of TSS
                keep = 0
                break
            if end >= v['tss'] and (start - v['tts']) < filter_tts:
                # start pos is within 20k of TTS
                keep = 0
                break
        elif v['strand'] == '-':
            if start <= v['tss'] and (v['tts'] - end) < filter_tts:
                keep = 0
                break
            if end >= v['tts'] and (start - v['tss']) < filter_tss:
                keep = 0
                break
    return keep

def sort_bed_like_file(fn):
    status = run_shell(f'sort -k 1,1 -k 2,2n {fn} > {fn}.sorttmp;mv {fn}.sorttmp {fn}')
    if status:
        logger.error(f'fail to sort file {fn}')
    return status

def change_enhancer(pwout, fn_count_enhancer, factors_d, n_ctrl, n_case, sam_ctrl, sam_case, flag=1):
    # flag = 1, means that the normalization factors are passed in, no need to calculate there
    condition = 2 if n_case > 0 else 1
    if flag != 1:
        logger.debug(f'nf passed in, but flag set to be {flag}, will ignore, and set flag as 1')
    
    data = pd.read_csv(fn_count_enhancer, sep='\t')
    data['Enhancer_ID'] = data['Enhancer_ID'].astype(str)
    fn_norm = f'{pwout}/eRNA/normalized_count_enhancer.txt'
    
    # Enhancer_ID	chr	start	end	4104-AW-1_sorted_rRNArm-F4q10	4104-AW-3_sorted_rRNArm-F4q10	4104-AW-2_sorted_rRNArm-F4q10	4104-AW-4_sorted_rRNArm-F4q10
    # 1	chr1	100080885	100091297	241	329	249	282
    # 2	chr1	103618441	103638394	3140	4570	3455	3746

    if condition == 1:
        # no case samples
        if n_ctrl == 1:
            try:
                force_symlink(fn_count_enhancer, fn_norm)
            except:
                logger.warning(f'fail to create symlink from {fn_count_enhancer} to {fn_norm}')
        else:
            pass # save norm data at the end
    else:
        # with case, will also calculate the Enhancer_change
        fn_change = f'{pwout}/eRNA/Enhancer_change.txt'
        
        if n_case + n_ctrl == 2:
            data_change = data.iloc[:, :4].copy()
            # each condition only have a single sample
            sam_ctrl, sam_case = sam_ctrl[0], sam_case[0]
            fc = ((data[sam_case] + 1) * factors_d[sam_case]) / ((data[sam_ctrl] + 1) * factors_d[sam_ctrl])
            data_change['log2fc'] = np.log2(fc)
            data_change.to_csv(fn_change, sep='\t', index=False, na_rep='NA')
        else:
            # run deseq2
            n_gene_cols = 4

            fn_metadata = f'{pwout}/intermediate/design_table_for_deseq2.txt'
            metadata = pd.read_csv(fn_metadata, sep='\t', index_col=0)

            ref_level = 'control'
            test_level = 'case'
            
            size_factors_in = 1 / np.array([factors_d[sam_lb] for sam_lb in sam_ctrl + sam_case])
            with HiddenPrints():
                res_df, size_factors = run_deseq2(n_gene_cols, data, metadata, ref_level, 
                size_factors_in=size_factors_in, test_level=test_level)
                res_df.to_csv(fn_change, sep='\t', index=False, na_rep='NA')

    for sam_lb, nf in factors_d.items():
        data[sam_lb] = data[sam_lb] * nf
    data.to_csv(fn_norm, sep='\t', index=False, na_rep='NA')



class Analysis:
    def __init__(self, args, is_long_eRNA=False, skip_get_mapped_reads=False, raw_input=True):
        
        self.bin_size = 200 # the bin size for building the pre-count dict of each bed file. compared bin size of 10, 20, 50, 100 and 200 and 500.  200 is the best, maybe because when getting the gb count, the size fit the distribution
        
        # folders
        self.bin_dir = os.path.dirname(os.path.realpath(__file__))
        self.pwout = args.pwout or os.getcwd()
        self.pwout_raw = args.pwout_raw
        self.pw_bed = f'{args.pwout_raw}/bed'
        self.inter_dir = os.path.join(self.pwout, 'intermediate')
        self.known_gene_dir = os.path.join(self.pwout, 'known_gene')
        self.longerna = is_long_eRNA # toggle if this is long eRNA
        self.skip_get_mapped_reads = skip_get_mapped_reads
        respect_sorted = vars(args).get('sorted', False)

        if in_docker:
            home = args.pwout_raw
        else:
            home = os.path.expanduser('~')
        self.home = home


        for lb, d in zip(['Output', 'Intermediate', 'Known Genes'], [self.pwout, self.inter_dir, self.known_gene_dir]):
            if not os.path.exists(d):
                logger.debug(f'Making {lb} directory: {d}')
                os.makedirs(d)

        self.status = 0

        self.organism = args.organism
        if self.organism.lower() not in {'hg19', 'hg38', 'mm10', 'dm3', 'dm6', 'ce10', 'danrer10', 'mm39'}:
            logger.warning(f"Non-built-in organism specified: {self.organism}, should comes with a valid GTF file and fasta file")

        # input files
        in1 = process_input(self.pwout_raw, args.in1, respect_sorted=respect_sorted) if raw_input else args.in1 # return = fn_lb, fn_bed
        if in1 is None:
            logger.error("Invalid input files provided for condition1")
            self.status = 1
        in2 = (process_input(self.pwout_raw, args.in2, respect_sorted=respect_sorted) or []) if raw_input else args.in2
        self.control_bed = in1
        self.case_bed = in2
        
        if self.longerna:
            self.n_gene_cols = 1
            self.gene_cols = ['long-eRNA']
            self.longerna_flag_str = '-longeRNA'
            self.longerna_prefix = 'longeRNA-'
        else:
            self.n_gene_cols = 2
            self.gene_cols = ['Transcript', 'Gene']
            self.longerna_flag_str = ''
            self.longerna_prefix = ''
        # ref files
        fa_in = vars(args).get('fa')
        if fa_in and os.path.exists(fa_in):
            if not fa_in.endswith('.fa'):
                logger.error(f"Invalid fasta file provided, file extension should be .fa in plain text. input = {fa_in}")
                self.status = 1
        ref_fls = get_ref(home, self.organism, fn_gtf=args.gtf, fn_fa=fa_in)
        if ref_fls is None:
            logger.error("Error encountered while retrieving reference files")
            self.status = 1
            return
        self.fa = ref_fls['fa']
        self.gtf = ref_fls['gtf']
        if not os.path.exists(self.gtf):
            logger.error(f"GTF file not available: {self.gtf}")
            self.status = 1
        
        # verfify the normalization factors
        factors = [
            [args.f1, self.control_bed, 'condition1'],
            [args.f2, self.case_bed, 'condition2']
        ]
        
        if self.case_bed:
            if bool(args.f1) != bool(args.f2):
                logger.error("Please provide normalization factors for both or none of the conditions")
                self.status = 1
        for factor, bed, lb in factors:
            n_bed = len(bed) if bed else 0
            if factor:
                if len(factor) != n_bed:
                    logger.error(f"Number of normalization factors provided for {lb} does not match the number of input files")
                    self.status = 1
                for n, f in enumerate(factor):
                    if f == 0:
                        logger.error(f"{lb} - factore {n + 1} / {len(factor)}: normalize factor could not be zero!")
                        self.status = 1
        
        if self.status:
            return
        
        self.ref = ref_fls
        self.input_fls = self.control_bed + (self.case_bed or []) # fn_lb, fn_bed
        
        fn_count_pp_gb = f'{self.inter_dir}/{self.longerna_prefix}count_pp_gb.txt'
        self.out_fls = {
            'bed_peaks': {},
            'fdr': {},
            'count_pp_gb': fn_count_pp_gb
        }
        
        

        
        fn_lb_uniq = set()
        fn_lb_dup = {}

        for fn_lb, fn_bed in self.input_fls:
            # tssCount, tssLength, tssRatio, tssSummit_pos, genebodyCount, genebodyLength, genebodyRatio, tss_vs_genebody_ratio
            if fn_lb not in fn_lb_uniq:
                fn_lb_uniq.add(fn_lb)
            else:
                fn_lb_dup.setdefault(fn_lb, []).append(fn_bed)
                continue
            header_bed_peaks = self.gene_cols + ['ppc', 'ppm', 'ppd', 'pps', 'gbc', 'gbm', 'gbd', 'pauseIndex']
            fn_bed_peaks = os.path.join(self.inter_dir, fn_lb + self.longerna_flag_str + '_raw.txt')

            if not (self.skip_get_mapped_reads and os.path.exists(fn_count_pp_gb)):
                fh_bed_peaks = open(fn_bed_peaks, 'w')
                fh_bed_peaks.write('\t'.join(header_bed_peaks) + '\n')
            else:
                fh_bed_peaks = open(fn_bed_peaks, 'r')
            self.out_fls['bed_peaks'][fn_lb] = {'fn': fn_bed_peaks, 'fh': fh_bed_peaks, 'header': header_bed_peaks}
            
            header_fdr = header_bed_peaks + ['pvalue', 'FDR']
            fn_fdr = os.path.join(self.inter_dir, fn_lb + self.longerna_flag_str + '_FDR.txt')
            self.out_fls['fdr'][fn_lb] = {'fn': fn_fdr, 'header': header_fdr}
        
        
        if fn_lb_dup:
            tmp = json.dumps(fn_lb_dup, indent=4)
            logger.error(f'Duplicate file labels found in input file list, please make the file name unique:\n{tmp}')
            sys.exit(1)
        
        # config
        self.config = {
            'mapped_sites_only': 0, # 1: only mapped sites are used to calculate the mappable sites, 0: all ATCG bases are used regardless of whether they are covered by reads. in the perl code, it is always set to 0
            'pro_up': args.up, # define the upstream of TSS as promoter (bp, default: 500)
            'pro_down': args.down, # define the downstream of TSS as promoter (bp, default: 500)
            'gb_start': args.gb, # define the start of gene body density calculation (bp, default: 1000)
            'min_gene_len': args.min_gene_len, # define the minimum length of a gene to perform the analysis (bp, default: 1000). Genes with length less than this will be ignored
            'window_size': args.window, # define the window size (bp, default: 50)
            'step': args.step, # define the step size (bp, default: 5)
            'normalize_factor1': args.f1, # normalization factors for samples of condition1
            'normalize_factor2': args.f2, # normalization factors for samples of condition2
            'tts_padding': args.tts_padding
        }
        

def get_gtf_intervals(gtf):
    intervals = {}
    for ts, v in gtf.items():
        chr_, start, end, strand, gene_name = v['chr'], v['start'], v['end'], v['strand'], v['gene_name']
        intervals.setdefault(chr_, []).append((start, end, gene_name, ts))
    for k, v in intervals.items():
        intervals[k] = sorted(v, key=lambda x: x[0])
    return intervals


def get_closest_downstream_gene(gtf, intervals, cumu_max, ts):
    # intervals = get_gtf_intervals(gtf)
    v = gtf[ts]
    chr_, start, end, strand, gene_name = v['chr'], v['start'], v['end'], v['strand'], v['gene_name']
    intervals_chr = intervals[chr_]
    cumu_max_chr = cumu_max[chr_]
    n_ts = len(intervals_chr)
    idx_left = bisect.bisect_left(cumu_max_chr, start)
    if strand == '+':
        for i in intervals_chr[idx_left:]:
            if i[2] != gene_name and i[1] > start:
                distance = i[0] - end
                return [*i, distance]
    else:
        idx = bisect.bisect_right(intervals_chr, (end,))
        if idx >= n_ts:
            idx = n_ts - 1
        min_distance = 1e9
        min_i = None
        n_cumu = len(cumu_max_chr)
        for i in intervals_chr[idx:0:-1]:
            try:
                i[2] != gene_name  and cumu_max_chr[idx] == i[1] and i[1] <= start
            except:
                print(f'ERROR, {ts} - {v} - n_cumu  = {n_cumu}, idx = {idx}, n_intervals = {len(intervals_chr)} - {i}')
            if i[2] != gene_name  and cumu_max_chr[idx] == i[1] and i[1] <= start:
                distance = start - i[1]
                if distance < min_distance:
                    min_distance = distance
                    min_i = i
                break
            if i[2] != gene_name and i[0] < end:
                distance = start - i[1]
                if distance < min_distance:
                    min_distance = distance
                    min_i = i
        if min_i is None:
            return [None, None, None, None, 1e-6]
        return [*min_i, min_distance]



def get_no_overlap_region_v2(gtf, func_get_target_region):
    # this is the best solution, only takes 120us, compared to 4.7s using dummy nested iteration
    # improved 40 times
    # func_get_target_region, is the function that get the target regioin from the gene_info
    #. e.g. get the downstream 50kb from TTS, 
    # input  = gene_info
    # return should be chr_, gn, s1, e1
    no_overlap = set()
    intervals = get_gtf_intervals(gtf)
    def get_cumu_max(intervals):
        # get the max end site position
        cumu_max = 0
        res = []
        for i in intervals:
            if i[1] > cumu_max:
                cumu_max = i[1]
            res.append(cumu_max)
        return res
    cumu_max = {k: get_cumu_max(v) for k, v in intervals.items()}
    for ts, v in gtf.items():
        chr_, gn, s1, e1 = func_get_target_region(v)
        overlap_found = False
        start_interval_idx = bisect.bisect_left(cumu_max[chr_], s1)
        start_interval_idx = bisect.bisect_left(cumu_max[chr_], cumu_max[chr_][start_interval_idx -1 if start_interval_idx > 0 else 0])
        for i in intervals[chr_][start_interval_idx:]:
            if i[2] == gn:
                continue
            if i[0] > e1:
                break
            if s1 <= i[1]:
                overlap_found = True
                break
        if not overlap_found:
            no_overlap.add(ts)
    return no_overlap
    

def process_bed_files(analysis, fls, gtf_info, gtf_info_raw, fa_idx, fh_fa, reuse_pre_count=False, save_tts_count=True, islongerna=False, downstream_no_overlap_length=50000):
    invalid_chr_transcript = 0
    bin_size = analysis.bin_size
    # ['ppc', 'ppm', 'ppd', 'pps', 'gbc', 'gbm', 'gbd', 'pauseIndex']
    fail_to_retrieve_seq = 0
    pwout = analysis.pwout
    window_size, step_size = analysis.config['window_size'], analysis.config['step']
    pw_bed = analysis.pw_bed
    n = 0
    prev = [time.time(), 0] # time, get_mapped_reads_count count
    section_size = 1000 # print log every 1000 loops

    # add more values to gtf_info
    s_before_loop = time.time()
    s = s_before_loop
    ts_list = list(gtf_info)

    pp_str = {ts: [] for ts in ts_list}
    gb_str = {ts: [] for ts in ts_list}
    tts_str = {ts: [str(gtf_info[ts][_]) for _ in ['chr', 'strand', 'tts']] for ts in ts_list}
    tts_down_str = {}
    downstream_no_overlap_length_str = f'{downstream_no_overlap_length/1000:.0f}k'
    
    # reformat the gtf_info by chrom
    gtf_info_new = {}  # k1 = chr, k2 = ts
    
    invlude_all_ts = True
    for ts, v in gtf_info.items():
        chr_ = v['chr']
        gtf_info_new.setdefault(chr_, {})[ts] = v
    if islongerna:
        ts_without_overlap = set()
    else:
        def get_down_tts_from_gene_info(gene_info):
            start, end, strand, chr_, gn = [gene_info[_] for _ in ['start', 'end', 'strand', 'chr', 'gene_name']]
            s1, e1 = (start, end + downstream_no_overlap_length) if strand == '+' else (start - downstream_no_overlap_length, end)
            return chr_, gn, s1, e1
        if not invlude_all_ts:
            ts_without_overlap = get_no_overlap_region_v2(gtf_info_raw, func_get_target_region=get_down_tts_from_gene_info)
            logger.debug(f'transcript without overlap from TSS to {downstream_no_overlap_length_str} downstream of TTS = {len(ts_without_overlap)}')
        else:
            ts_without_overlap = set()
    
    # logger.warning(f'modify here, now will get the {downstream_no_overlap_length_str} downstream count for all transcripts')
    
    
    # seq_pool = {} # key = transcript_id
    # count_pool = {}
    # prev_peak_pool = {}
    # genes_with_N = set()
    n_ts_init = len(ts_list)
    err = 0
    for fn_lb, fn_bed in fls:
        logger.info(f'processing {fn_lb}')
        fh_bed_peaks = analysis.out_fls['bed_peaks'][fn_lb]['fh']

        pre_count_flist = pre_count_for_bed(fn_lb, fn_bed, pw_bed, bin_size, reuse=reuse_pre_count)
        
        prev_peak = {}
        # exclude the ts that the chr not in the bed files
        valid_chr = set(pre_count_flist)
        chr_excluded = {}
        ts_excluded = 0
        
        for ts, v in gtf_info.items():
            ts_chr = v['chr']
            if ts_chr not in valid_chr:
                ts_excluded += 1
                gtf_info_new[ts_chr].pop(ts)
                if ts in pp_str:
                    del pp_str[ts]
                if ts in gb_str:
                    del gb_str[ts]
                if ts in tts_str:
                    del tts_str[ts]
                chr_excluded.setdefault(ts_chr, 0)
                chr_excluded[ts_chr] += 1
        if ts_excluded > 0:
            logger.info(f'{ts_excluded} transcripts in GTF file excluded due to chromosome not exist in bed file')
            logger.debug(f'excluded chr = {chr_excluded}')
        logger.debug(f'init ts list = {n_ts_init}, current = {n_ts_init - ts_excluded}, drop = {ts_excluded}')

        n_current = n_ts_init - ts_excluded
        if n_current == 0:
            size = os.path.getsize(fn_bed)
            if size < 100:
                logger.error(f'the input bed file is empty, skip')
                err = 1
            else:
                logger.error(f'no transcript left after matching chrom between gtf and bed file, maybe because of different chrom name pattern, skip')
                err = 1

        for chr_, ts_list_chr in gtf_info_new.items():
            fn_count_per_base, fn_count_bin = pre_count_flist[chr_]
            with open(fn_count_per_base, 'rb') as f:
                count_per_base = pickle.load(f)
            with open(fn_count_bin, 'rb') as f:
                count_bin = pickle.load(f)
            bin_size = count_bin['bin_size']
            
            for transcript_id in ts_list_chr:
                gene_info = gtf_info[transcript_id]
                chr_, strand, gene_raw_s, pp_start, pp_end, gb_start, gb_end, strand_idx, gb_len_mappable, gene_seq, tts_start, tts_end, tts_down_region_s, tts_down_region_e = [gene_info[_] for _ in ['chr', 'strand', 'start', 'pp_start', 'pp_end', 'gb_start', 'gb_end', 'strand_idx', 'gb_len_mappable', 'gene_seq', 'tts_start', 'tts_end', 'tts_down_region_s', 'tts_down_region_e']]

                n += 1
                if n % section_size == 0:
                    now = time.time()
                    prev_time, prev_count = prev
                    time_gap = now - prev_time
                    speed = time_gap * 1000  / (n - prev_count)
                    logger.debug(f'{n/1000}k - get_mapped_reads_count count, time_gap={time_gap:.2}s, speed={speed:.3f}ms')
                    prev = [now, n]

                # pp_res: {'ppc': prev[0], 'ppd': ppd, 'mappable_sites': mappable_sites, 'summit_pos': summit_pos_str, 'summit_count': summit_count}
                s = time.time()
                gbc, gbd, pp_res, tts_ct = get_peak(count_per_base, count_bin, chr_, strand, gene_raw_s, strand_idx, pp_start, pp_end, gb_start, gb_end, tts_start, tts_end, gb_len_mappable, gene_seq,  window_size, step_size, bin_size, prev_peak)
                # get_peak_method2(count_per_base, count_bin, chr_, strand_idx, s, e, bin_size)
                # get the count for last exon and the downstream region
                if invlude_all_ts or transcript_id in ts_without_overlap:
                    last_exon_s, last_exon_e = gene_info['last_exon']
                    last_exon_ct = get_peak_method2(count_per_base, count_bin, chr_, strand_idx, last_exon_s, last_exon_e, bin_size)
                    tts_down_ct = get_peak_method2(count_per_base, count_bin, chr_, strand_idx, tts_down_region_s, tts_down_region_e, bin_size)
                    last_exon_len = last_exon_e - last_exon_s + 1
                    if transcript_id not in tts_down_str:
                        tts_down_str[transcript_id] = [transcript_id, gene_info['gene_name'], chr_, strand, str(last_exon_len), str(last_exon_s), str(last_exon_e)]
                    tts_down_str[transcript_id] += [str(last_exon_ct), str(tts_down_ct), f'{tts_down_ct / last_exon_ct:.4f}' if last_exon_ct > 0 else 'NA']

                pp_str[transcript_id].append(str(pp_res['ppc']))
                gb_str[transcript_id] += [str(gbc), str(gbd)]
                tts_str[transcript_id].append(str(tts_ct))
                
                # write to file
                # row = ['Transcript', 'Gene', 'tssCount', 'tssLength', 'tssRatio', 'tssSummit_pos', 'genebodyCount', 'genebodyLength', 'genebodyRatio', 'tss_vs_genebody_ratio']
                row = [transcript_id]
                if not analysis.longerna:
                    row.append(gene_info['gene_name'])
                row += [pp_res['ppc'], pp_res['mappable_sites'], pp_res['ppd'], pp_res['summit_pos'], gbc, gb_len_mappable, gbd]
                pro_vs_pb = round(pp_res['ppd'] / gbd, 4) if gbd else 'NA'
                row.append(pro_vs_pb)
                
                print('\t'.join(map(str, row)), file=fh_bed_peaks)
            # # release memory
            del count_per_base
            del count_bin 
            gc.collect()
    
    if err:
        sys.exit(1)
    if fail_to_retrieve_seq:
        logger.warning(f'failed to retrieve sequence for {fail_to_retrieve_seq} genes')

    if invalid_chr_transcript:
        logger.warning(f"Number of genes with invalid chromosome = {invalid_chr_transcript}")

    # save tts down 50k results
    fn_tts_down = f'{pwout}/intermediate/tts_down_{downstream_no_overlap_length_str}.txt'
    logger.info(f'number of transcripts {downstream_no_overlap_length_str} downstream region to file = {len(tts_down_str)}, ts_without_overlap_{downstream_no_overlap_length_str}={len(ts_without_overlap)}')
    with open(fn_tts_down, 'w') as o:
        header = ['Transcript', 'Gene', 'chr', 'strand', 'last_exon_len', 'last_exon_s', 'last_exon_e']
        for fn_lb, _ in fls:
            header += [f'last_exon_{fn_lb}', f'tts_down_{downstream_no_overlap_length_str}_{fn_lb}', f'ratio_{fn_lb}']
        print('\t'.join(header), file=o)
        for ts, v in tts_down_str.items():
            print('\t'.join(v), file=o)

    # save the tts_count
    if save_tts_count:
        fn_tts_count = f'{pwout}/intermediate/count_tts.raw.txt'
        fn_tts_count_filtered = f'{pwout}/intermediate/count_tts.filtered.txt'
        tts_file_header = ['Transcript', 'Gene', 'chr', 'strand', 'TTS'] + [_[0] for _ in fls]

        # filter count_tts, exclude the trascripts which have other TSS within 2k downstream region
        # {'chr': '1',
        # 'strand': '+',
        # 'gene_name': 'SGIP1',
        # 'start': 66999639,
        # 'end': 67216822}
        # first get all TSS as a sorted list
        tss_list = {}  # key is chr
        gene_regions = {}
        for gene_info in gtf_info.values():
            chr_ = gene_info['chr']
            tss_list.setdefault(chr_, set()).add(gene_info['tss'])
            gene_regions.setdefault(chr_, []).append([gene_info['start'], gene_info['end'], gene_info['gene_name']])
        tss_list = {k: sorted(v) for k, v in tss_list.items()}
        gene_regions = {k: sorted(v, key=lambda x: x[0]) for k, v in gene_regions.items()}
        gene_regions_start = {k: [_[0] for _ in v] for k, v in gene_regions.items()}
        gene_regions_max_len = {k: max([_[1] - _[0] + 1 for _ in v]) for k, v in gene_regions.items()}

        tts_dowstream_no_tss_region_len = 2000  # exclude the trascripts which have other TSS within 2k downstream region
        n_excluded = 0
        # also, there shouldn't be any overlap between TTS and any gene region
        # during exluding the TTS with any other gene regions, 
        # for hg19, 5255 transcripts are excluded (15.6%), 28519 are kept
        
        with open(fn_tts_count, 'w') as o, open(fn_tts_count_filtered, 'w') as o1:
            print('\t'.join(tts_file_header), file=o)
            print('\t'.join(tts_file_header), file=o1)
            for ts in ts_list:
                gene_info = gtf_info[ts]
                gn = gene_info['gene_name']
                if ts in tts_str:
                    line = '\t'.join([ts, gn] + tts_str[ts])
                    print(line, file=o)
                    
                    
                    tts = gene_info['tts']
                    strand = gene_info['strand']
                    chr_ = gene_info['chr']
                    tss_list_chr = tss_list[chr_]
                    # try:
                    tts_region_s, tts_region_e = [tts, tts + tts_dowstream_no_tss_region_len] if strand == '+' else [tts - tts_dowstream_no_tss_region_len, tts]
                    left_idx =  bisect.bisect_left(tss_list_chr, tts_region_s)
                    right_idx = bisect.bisect_right(tss_list_chr, tts_region_e)
                    # except:
                    #     e = traceback.format_exc()
                    #     logger.error(e)
                    #     logger.info([tts, tts_dowstream_no_tss_region_len])
                    #     logger.info(f'chr = {chr_}, tts_region_s={tts_region_s}, tts_region_e= {tts_region_e}, tss_list_chr len = {len(tss_list_chr)}, first 20 = {tss_list_chr[:20]}' )
                    #     sys.exit(1)
                        
                    if left_idx < right_idx:
                        n_excluded += 1
                        continue
                    
                    # check if there are any overlap with other genes
                    idx_tts_in_start_pos = bisect.bisect_right(gene_regions_start[chr_], tts)
                    # iterate over the left section
                    excluded = False
                    gene_regions_max_len_chr = gene_regions_max_len[chr_]
                    
                    # here must use idx_tts_in_start_pos, ensure all regions iterate with start pos <= tts
                    for iregion in gene_regions[chr_][:idx_tts_in_start_pos][::-1]:
                        if iregion[1] >= tts:
                            if gn == iregion[2]:
                                # same gene, skip
                                continue
                            n_excluded += 1
                            excluded = 1
                            break
                        if tts - iregion[0] > gene_regions_max_len_chr:
                            # no overlap
                            break
                    if excluded:
                        continue
                    print(line, file=o1)
        
        logger.debug(f'TTS count file saved: {fn_tts_count}, filtered = {fn_tts_count_filtered}, transcript excluded during filtering due to having TSS within 2k downstream of TTS region = {n_excluded}')
    
    return pp_str, gb_str

def get_new_args(args, update_dict):
    args_dict = vars(args)
    args_dict.update(update_dict)
    return SimpleNamespace(**args_dict)


def parse_design_table(args):
    pwout_raw = os.path.realpath(args.pwout)
    args.pwout = pwout_raw
    args.pwout_raw = pwout_raw
    args.pw_bed = f'{pwout_raw}/bed'
    if not os.path.exists(args.pw_bed):
        os.makedirs(args.pw_bed, exist_ok=True)
    
    if args.gtf:
        args.gtf = os.path.realpath(args.gtf)

    if args.in1 is not None:
        return [['', args]]
    elif args.design_table is not None:
        group_info = {}  # key = group name, v = file list
        comparison = []
        groups_in_comparison = set()
        group_info['null'] = []
        new_arg_list = []
        err = 0
        with open(args.design_table) as f:
            batch_map = {}
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                if line.startswith('@@'):
                    tmp = re.split(r'\s+', line[2:].strip())
                    if len(tmp) != 2:
                        logger.error(f"Invalid comparison definition: {line}, should have 2 columns, col1 = case group, col2 = control group")
                        sys.exit(1)
                    comparison.append(tmp)
                    groups_in_comparison |= set(tmp)
                else:
                    tmp = re.split(r'\s+', line.strip())
                    if len(tmp) not in {2, 3}:
                        logger.error(f"Invalid line: {line}, should have 2 or 3 columns, col1 = file path, col2 = group name, col3 = optional, the batch of this sample")
                        err = 1
                        continue
                    fn = os.path.realpath(tmp[0])
                    if not os.path.exists(fn):
                        logger.error(f'input file not exist: {tmp[0]}')
                        err = 1
                        continue
                    group_info.setdefault(tmp[1], []).append(fn)
                    if len(tmp) == 3:
                        batch_map[fn] = tmp[2] or 'empty'
            total_fls = sum([len(_) for _ in group_info.values()])
            if batch_map and len(batch_map) != total_fls:
                logger.error(f'Not all samples have the batch information (3rd column), total files = {total_fls}, with batch info = {len(batch_map)}, without batch info = {total_fls - len(batch_map)}')
                err = 1
            if err:
                logger.error('quit now.')
                sys.exit(1)
            group_not_defined = groups_in_comparison - set(group_info)
            if group_not_defined:
                logger.error(f"Groups in comparison not defined in the sample section: {sorted(group_not_defined)}, groups defined in sample = {sorted(group_info)} quit now")
                sys.exit(1)
            
            # if no comparison defined, process each group separately
            if len(comparison) == 0:
                logger.error(f'No comparison defined in {args.design_table}, please use @@ to define the comparison .If need to process a group separately, use @@null\\tgroup_name, quit now.')
                sys.exit(1)
            else:
                for case, control in comparison:
                    update_dict = {}
                    update_dict['in1'] = group_info[control]
                    update_dict['in2'] = group_info[case] or None
                    if batch_map:
                        all_sams = [_ for _ in update_dict['in1']]
                        if update_dict['in2']:
                            all_sams += [_ for _ in update_dict['in2']]
                        update_dict['batches'] = [batch_map[_] for _ in all_sams]
                    else:
                        update_dict['batches'] = None
                    comp_str = control if case == 'null' else f'{case}_vs_{control}'
                    pwout = os.path.join(pwout_raw, comp_str)
                    update_dict['pwout'] = pwout
                    
                    args_new = get_new_args(args, update_dict)
                    new_arg_list.append([comp_str, args_new])
        return new_arg_list
    else:
        logger.error(f"No input files provided, either use -in1 / in2 or provide a design table, input arguments = {vars(args)}")
        sys.exit(1)


def check_full_rank(condition, batches):
    """
    if the batch is fully correlated with condition, then, there is no need to add the batch column
    """
    condition = [str(_) for _ in condition]
    batches = [str(_) for _ in batches]
    df = pd.DataFrame({'condition': condition, 'batch': batches})
    # use drop_first=True
    design_matrix = pd.get_dummies(df, drop_first=True)
    design_matrix = design_matrix.astype(float)
    rank = np.linalg.matrix_rank(design_matrix.values)
    return rank == design_matrix.shape[1]



def build_design_table(pwout, in1, in2, batches):
    fn_metadata = f'{pwout}/intermediate/design_table_for_deseq2.txt'
    logger.debug(f'building design table for deseq2: {fn_metadata}')
    in2 = in2 or []
    fls = in1 + in2
    fls_lb = [get_lb(_) for _ in fls]
    rep1, rep2 = len(in1), len(in2)
    condition = ['control'] + ['case'] * (rep1 - 1) if rep2 == 0 else ['control'] * rep1 + ['case'] * rep2
    if batches and len(batches) != len(fls):
        logger.error(f'batch number ({len(batches)}) not match with total sample number ({len(fls)}! exit')
        sys.exit(1)
        
    # check for fully ranked
    if batches:
        is_full_rank = check_full_rank(condition, batches)
        if not is_full_rank:
            logger.warning(f'Batch information provided, however, the model matrix is not full rank, condition column is fully correlated with batch column, will remove the batch column.')
            batches = None
            
        
    with open(fn_metadata, 'w') as o:
        if batches:
            print(f'\tcondition\tbatch', file=o)
            for lb, icond, ibatch in zip(fls_lb, condition, batches):
                print(f'{lb}\t{icond}\t{ibatch}', file=o)
        else:
            print(f'\tcondition', file=o)
            for lb, icond in zip(fls_lb, condition):
                print(f'{lb}\t{icond}', file=o)


def pause_longeRNA_main(args):
    """
    args: an object with attributes equiv to the argparse object, must specify 
    """

    # check if the attributes are present
    defined_attrs = vars(args)
    required_attrs = {'in1', 'pwout', 'organism', 'pw_bed', 'gtf'}
    # f1, f2 = normalization factors
    
    optional_attrs = {
        'in2': None,
        'fa': None,
        'up': 500,
        'down': 500,
        'gb': 1000,
        'min_gene_len': 1000,
        'window': 50,
        'step': 5,
        'ignore': False,
        'skip_get_mapped_reads': False,
    }
    pwout = args.pwout
    pw_bed = args.pw_bed
    args.tts_padding = 2000
    missing = required_attrs - set(defined_attrs)
    exist = required_attrs & set(defined_attrs)
    for attr in exist:
        if getattr(args, attr) is None:
            missing.add(attr)
    if missing:
        logger.error(f"Missing required attributes: {sorted(missing)}")
        return 1
    for attr, default in optional_attrs.items():
        if attr not in defined_attrs:
            setattr(args, attr, default)

    if in_docker:
        home = args.pwout_raw
    else:
        home = os.path.expanduser('~')
    

    # check dependencies
    skip_get_mapped_reads = args.skip_get_mapped_reads
    dependency_status = check_dependency()
    if dependency_status:
        logger.error("Dependency check failed")
        sys.exit(1)

    analysis = Analysis(args, is_long_eRNA=True, raw_input=False, skip_get_mapped_reads=skip_get_mapped_reads) # raw_input set as False, will not run the process_input function
    if analysis.status:
        logger.error("Exit now")
        sys.exit(1)
    
    logger.debug(f'running pause_longeRNA: {vars(args)}')

    rep1 = len(analysis.control_bed)
    rep2 = len(analysis.case_bed) if analysis.case_bed else 0
    window_size = analysis.config['window_size']
    
    # get the normalization factors
    factor_flag = 1
    fn_nf = f'{pwout}/intermediate/nf.txt'
    nf_d = {}
    try:
        with open(fn_nf) as f:
            f.readline() # header
            for i in f:
                lb, v = i[:-1].split('\t')
                nf_d[lb] = float(v)
        logger.debug(f'previous normalization factors loaded')
    except:
        e = traceback.format_exc()
        logger.error(f'fail to get normalization factor, please make sure pause_PROseq is executed')
        logger.error(e)
        sys.exit(1)
    
    factor1 = [nf_d[lb] for lb, _ in analysis.control_bed]
    factor2 = [nf_d[lb] for lb, _ in analysis.case_bed]
    
    logger.debug(f'nf = {nf_d}')

    pro_up, pro_down, gb_down_distance, min_gene_len, tts_padding = [analysis.config[_] for _ in ['pro_up', 'pro_down', 'gb_start', 'min_gene_len', 'tts_padding']]

    # process the GTF file
    fn_gtf = analysis.ref['gtf']
    logger.debug(f"Processing longeRNA gtf file: {fn_gtf}")
    
    gtf_info = {}
    with open(fn_gtf) as f:
        f.readline()  # header
        # chr     start   end     strand
        # chr1    27365605        27376193        -
        # chr1    99872187        99904465        -
        # chr1    103653695       103664771       -
        # {'chr': chr_, 'strand': strand, 'gene_name': gene_name, 'start': start, 'end': end}
        for i in f:
            chr_, s, e, strand = i[:-1].split('\t')
            ts = f'{chr_}:{s}-{e}:{strand}'
            s, e = int(s), int(e)
            tss, tts = (s, e) if strand == '+' else (e, s)
            gene_name = f'{chr}:{s}-{e}:{strand}'
            ires = {'chr': chr_, 'strand': strand, 'start': s, 'end': e, 'tss': tss, 'tts': tts, 'gene_name': gene_name}
            gtf_info[ts] = add_value_to_gtf(ires, pro_up, pro_down, gb_down_distance, tts_padding, islongerna=True) 
    
    # prepare fasta file
    fn_fa = analysis.ref['fa']
    fa_idx = build_idx_for_fa(home, fn_fa)
    fh_fa = open(fn_fa, 'r')

    fls = analysis.input_fls # element is [fn_lb, fn_bed]
    sam_order = [_[0] for _ in fls]
    n_gene_cols = analysis.n_gene_cols

    def close_fh():
        # close file handle
        for fn_lb, fn_bed in fls:
            analysis.out_fls['bed_peaks'][fn_lb]['fh'].close()
        fh_fa.close()

    fn_count_pp_gb = analysis.out_fls['count_pp_gb']
    if skip_get_mapped_reads and os.path.exists(fn_count_pp_gb):
        logger.warning(f'demo mode, skip getting pp_gb count, reuse prev res')
        close_fh()
    else:
        reuse_pre_count = True
        logger.info(f'Getting pp_gb count')
        pp_str, gb_str = process_bed_files(analysis, fls, gtf_info, gtf_info, fa_idx, fh_fa, reuse_pre_count=reuse_pre_count, save_tts_count=False, islongerna=True)
        close_fh()
    
        # count_pp_gb
        logger.info('Building count_pp_gb')
        prefix = 'longeRNA-' if analysis.longerna else ''
        header = analysis.gene_cols
        header_pp = []
        header_gb = []
        for fn_lb, fn_bed in fls:
            header_gb += [f'gbc_{fn_lb}', f'gbd_{fn_lb}']
            header_pp.append(f'ppc_{fn_lb}')
        header += header_pp + header_gb

        
        with open(fn_count_pp_gb, 'w') as o:
            print('\t'.join(header), file=o)
            for transcript_id in pp_str:
                # gene_info = gtf_info[transcript_id]
                row = [transcript_id]
                row += pp_str[transcript_id] + gb_str[transcript_id]
                print('\t'.join(row), file=o)


    # calculate FDR
    close_fh()
    logger.info('Calculating FDR')
    pause_index_str = {}
    for fn_lb, fn_bed in fls:
        fn_peak = analysis.out_fls['bed_peaks'][fn_lb]['fn']
        fn_fdr = analysis.out_fls['fdr'][fn_lb]['fn']
        pause_index_str = get_FDR_per_sample(fn_lb, fn_peak, fn_fdr,  pause_index_str)

    pause_index_str_new = {}
    for transcript_id, v1 in pause_index_str.items():
        v2 = []
        for isam in sam_order:
            v3 = v1.get(isam, ['NA', 'NA', 'NA'])
            v2 += v3
        pause_index_str_new[transcript_id] = v2
    pause_index_str = pause_index_str_new
    

    fno_prefix = 'longeRNA-' if analysis.longerna else ''

    logger.info('Change_pp_gb')
    change_pp_gb(n_gene_cols, fn_count_pp_gb, analysis.pwout, rep1, rep2, window_size, factor1=factor1, factor2=factor2, factor_flag=factor_flag, islongerna=True)
    
    logger.debug('Dump pindex.txt')
    header_extra = []
    for fn_lb, fn_bed in fls:
        header_extra += [f'{fn_lb}-pindex', f'{fn_lb}-pvalue', f'{fn_lb}-FDR']
    header = analysis.gene_cols + header_extra
    fno =  f'{pwout}/eRNA/{fno_prefix}pindex.txt'
    fno_pindex_change = f'{pwout}/eRNA/{fno_prefix}pindex_change.txt'
    with open(fno, 'w') as o:
        print('\t'.join(header), file=o)
        for transcript_id, data in pause_index_str.items():
            row = [transcript_id]
            if not analysis.longerna:
                row.append(gtf_info[transcript_id]['gene_name'])
            row += data
            print('\t'.join(map(str, row)), file=o)

    if rep2 > 0:
        # change_pindex
        logger.info('Running change_pindex')
        change_pindex(fno_prefix, n_gene_cols, fn_count_pp_gb, fno_pindex_change, rep1, rep2, window_size, factor1=factor1, factor2=factor2, factor_flag=factor_flag)


def get_alternative_isoform_across_conditions(fn, pwout, pw_bed, rep1, rep2, tts_padding=None):
    """
    get the alternative isoform across conditions
    fn = count_pp_gb.txt / count_tts.txt
    the normalized version
    """
    n_transcript = {}
    data = pd.read_csv(fn, sep='\t')
    n_orig = len(data)
    n_transcript['input'] = n_orig
    
    # with open(fn_active_gene) as f:
    #     active_ts_set = {_.split('\t')[0] for _ in f}
    # logger.info(f'active transcripts = {len(active_ts_set)}')
    
    # keep active transcripts only
    # !!!! this is wrong , because the active_genes already collapsed by the gene, only 1 transcript for each gene
    
    cols_orig = list(data.columns)
    if 'count_pp_gb' in fn:
        sam_list = [_[4:] for _ in cols_orig if _.startswith('ppc_')]
        data.columns = [_[4:] if _.startswith('ppc_') else _ for _ in cols_orig]
        out_prefix = 'TSS_'
        pos_in_use = 'TSS'
        data[pos_in_use] = 0
        data.loc[data.strand == '+', pos_in_use] = data.loc[data.strand == '+', 'start']
        data.loc[data.strand == '-', pos_in_use] = data.loc[data.strand == '-', 'end']
        min_sum = 50
        distance_thres = 1000 # previous is 500
    elif 'count_tts' in fn:
        out_prefix = 'TTS_'
        pos_in_use = 'TTS'
        min_sum = 80
        distance_thres = 2000 # previous is 1k
        tmp = list(data.columns)
        idx_tts = tmp.index('TTS')
        sam_list = tmp[idx_tts+1:]
        if tts_padding is None:
            logger.error(f'tts_padding is not passed in for TTS input: {fn}')
            return 1
    else:
        logger.error(f'Input file should be count_pp_gb or count_tts, input = {fn}')
        return 1

    if rep2 < 1:
        logger.warning(f'No case samples, skip finding alternative isoforms')
        return 1
    
    # round the value
    round_digit = 3
    data[sam_list] = data[sam_list].round(round_digit)
    # only keep the genes with multiple isoforms with different TSS
    data = data.groupby('Gene').filter(lambda x: x['Transcript'].nunique() > 1)
    n_transcript['keep_only_multiple_isoforms'] = len(data)
    data = data.groupby('Gene').filter(lambda x: x[pos_in_use].nunique() > 1) # 333 rows
    n_transcript[f'remove_all_transcripts_with_same_{pos_in_use}'] = len(data)

    sam1, sam2 = (sam_list[:rep1], sam_list[rep1:])
    
    data['ctrl_sum'] = data[sam1].sum(axis=1)
    data['case_sum'] = data[sam2].sum(axis=1)
    
    # exclude transcripts with ctrl_sum and case_sum are both < 5
    if pos_in_use == 'TSS':
        data = data.loc[(data['ctrl_sum'] >= min_sum) | (data['case_sum'] >= min_sum)]
        n_transcript[f'filter_by_count_{pos_in_use}'] = len(data)
    elif pos_in_use == 'TTS':
        # for TTS, use the total sum of gbc as the normalization factor
        # below is the previous change_pp_gb code
        # the sum is about 37%-50% of the total reads count
        fn_gbc_sum = f'{pwout}/intermediate/gbc_sum.json'
        # tts_region_len = tts_padding * 2 # both upstream and downstream
        tts_region_len = tts_padding # currently, only use the downstream region
        with open(fn_gbc_sum) as f:
            gbc_sum = json.load(f)
        norm_density = data[sam_list]/gbc_sum * 10_000_000 / tts_region_len  # divide by region length to get density
        norm_density_mean = norm_density.sum(axis=1)/len(sam_list)
        
        data = data.loc[norm_density_mean > 0.0004]
        n_transcript[f'filter_by_count_{pos_in_use}'] = len(data)
        

        # for sn, idx_gbc_tmp, idx_gbd_tmp in zip(range(n_sam), idx_gbc_combined, idx_gbd_combined):
        #     sum_col = data.iloc[:, idx_gbc_tmp].sum()
        #     data[f'tmp_gbd_pass_{sn+1}'] = data.iloc[:, idx_gbd_tmp] * 10_000_000 / sum_col
        # cols_gbd_pass = [f'tmp_gbd_pass_{sn+1}' for sn in range(n_sam)]
        # data['tmp_gbd_pass'] = data[cols_gbd_pass].apply(lambda x: x.sum()/n_sam > 0.004, axis=1)
        
        # get the gbc sum for each sample
    else:
        logger.error(f'unknown input for identifying isoform: {fn}')
        sys.exit(1)
    
    # keep only genes still have multiple isoforms
    data = data.groupby('Gene').filter(lambda x: x['Transcript'].nunique() > 1)
    data = data.groupby('Gene').filter(lambda x: x[pos_in_use].nunique() > 1) # 333 rows
    n_transcript[f'keep_only_genes_with_multiple_isoforms_after filtering_count'] = len(data)
    
    tmp = json.dumps(n_transcript, indent=3)
    logger.debug(f'g@{pos_in_use}\n{tmp}')

    
    
    n_rows = len(data)
    if n_rows < 2:
        logger.warning(f'No transcrits remain for alternative isoform identification')
        return 1

    data = data.sort_values(['chr', 'Gene', pos_in_use])
    data['ctrl_mean'] = (data['ctrl_sum'] / rep1).round(round_digit)
    data['case_mean'] = (data['case_sum'] / rep2).round(round_digit)
    
    def compare_transcripts(gn, ts1, ts2, ctrl_max, case_max):
        ts_id1, ts_id2 = ts1.Transcript, ts2.Transcript
        tss_pos_ts1, tss_pos_ts2 = [ts1[pos_in_use], ts2[pos_in_use]]
        ratio_ctrl = (ts1['ctrl_mean'] / ts2['ctrl_mean']).round(round_digit) if ts2['ctrl_mean'] > 0 else np.nan
        ratio_case = (ts1['case_mean'] / ts2['case_mean']).round(round_digit) if ts2['case_mean'] > 0 else np.nan
        ts1_mean_ctrl, ts1_mean_case = ts1['ctrl_mean'], ts1['case_mean']
        ts2_mean_ctrl, ts2_mean_case = ts2['ctrl_mean'], ts2['case_mean']
        
        isoform_switch_flag = ''
        main_ts = {'case': set(), 'ctrl': set()}
        for condition, max_v in zip(['ctrl', 'case'], [ctrl_max, case_max]):
            for ts_lb, its in [['ts1', ts1], ['ts2', ts2]]:
                # logger.info(its)
                # logger.info(f'{condition}_mean')
                # logger.info(its['ctrl_mean'])
                v = its[f'{condition}_mean']
                if v == max_v:
                    main_ts[condition].add(ts_lb)
        main_ts_ctrl, main_ts_case = main_ts['ctrl'], main_ts['case']
        if (main_ts_ctrl and main_ts_case) and len(main_ts_case & main_ts_ctrl) == 0:
            isoform_switch_flag = 'isoform_switched'
        tables = []
        
        # if case_sum or ctrl_sum are zero for both ts, will not be suitable for cmhtest
        if ts1['ctrl_sum'] + ts2['ctrl_sum'] == 0 or ts1['case_sum'] + ts2['case_sum'] == 0:
            pvalue = np.nan
            odds_ratio = 0
        else:
            for isam1 in sam1:
                for isam2 in sam2:
                    tables.append([[int(ts1[isam1]), int(ts1[isam2])], [int(ts2[isam1]), int(ts2[isam2])]])
            if len(tables) == 0:
                pvalue = odds_ratio = np.nan
            elif len(tables) == 1:
                vals = tables[0][0] + tables[0][1]
                pvalue = fisher.pvalue(*vals).two_tail
                odds_ratio = ratio_ctrl / ratio_case
            else:
                with warnings.catch_warnings(record=True) as w:
                    odds_ratio, pvalue = cmhtest(tables)
                    if w:
                        logger.debug(f'cmhtest warning: {w[0].message}')
                        logger.debug(f'gn = {gn}, ts1 = {ts_id1} ts2 = {ts_id2}, tables = {tables}')
        return [ts1['chr'], gn, ts_id1, ts_id2, tss_pos_ts1, tss_pos_ts2, ctrl_max, ts1_mean_ctrl, ts2_mean_ctrl, ratio_ctrl, case_max, ts1_mean_case, ts2_mean_case, ratio_case, odds_ratio, pvalue, ';'.join(sorted(main_ts_ctrl)), ';'.join(sorted(main_ts_case)), isoform_switch_flag]

    def parse_by_gene(gn, g):
        """
        input is a groupby obj, with a single gene
        each row in the output is a transcript comparison.
        e.g. if there are 4 transcripts, the output will have 4*3/2 = 6 combinations
        """
        ts_list = list(g['Transcript'])
        n_ts = len(ts_list)
        res = []
        ctrl_max = g.ctrl_mean.max()
        case_max = g.case_mean.max()
        for i in range(n_ts - 1):
            ts1 = g.iloc[i]
            pos_ts1 = ts1[pos_in_use]
            for j in range(i + 1, n_ts):
                ts2 = g.iloc[j]
                pos_ts2 = ts2[pos_in_use]

                # if the TSS / TTS of these 2 transcripts are too close, skip this combination
                if abs(pos_ts1 - pos_ts2) < distance_thres:
                    continue
                res.append(compare_transcripts(gn, ts1, ts2, ctrl_max, case_max))
        return pd.DataFrame(res, columns=['chr', 'Gene', 'Transcript1', 'Transcript2', f'{out_prefix}transcript1', f'{out_prefix}transcript2', 'Ctrl_max_count', 'Ctrl_ts1_count', 'Ctrl_ts2_count', 'Ratio_ts1_vs_ts2_in_Ctrl', 'Case_max_count', 'Case_ts1_count', 'Case_ts2_count', 'Ratio_ts1_vs_ts2_in_Case', 'odds_ratio', 'pvalue', 'main_isoform_ctrl', 'main_isoform_case', 'isoform_switched'])
    
    
    df_isoform = None
    empty_df = []
    for gn, g in data.groupby('Gene'):
        tmp = parse_by_gene(gn, g)
        if len(tmp) == 0:
            empty_df.append(gn)
            continue
        if df_isoform is None:
            df_isoform = tmp
        else:
            df_isoform = pd.concat([df_isoform, tmp])
        
    if len(empty_df) > 0:
        logger.debug(f'gn without transcripts for comparison due to {pos_in_use} too close, n = {len(empty_df)}, first 10 genes = {empty_df[:10]}')
    
    df_isoform = add_FDR_col(df_isoform, 'pvalue')
    df_isoform.sort_values(['isoform_switched', 'pvalue', 'chr', 'Gene'], ascending=[False, True, True, True], inplace=True)
    sig = df_isoform.loc[df_isoform.pvalue < 0.05]
    
    fn_isoform = f'{pwout}/intermediate/{out_prefix}alternative_isoforms_between_conditions.tsv'
    fn_sig = f'{pwout}/known_gene/{out_prefix}alternative_isoforms_between_conditions.sig.tsv'
    df_isoform.to_csv(fn_isoform, index=False, na_rep='NA', sep='\t')
    sig.to_csv(fn_sig, index=False, na_rep='NA', sep='\t')
    
    # export excel format, sometimes, the package dependency for excel is not installed, so use try except
    fn_isoform_xlsx = fn_isoform.replace('.tsv', '.xlsx')
    fn_sig_xlsx = fn_sig.replace('.tsv', '.xlsx')
    try:
        df_isoform.to_excel(fn_isoform_xlsx, index=False, na_rep='NA')
        sig.to_excel(fn_sig_xlsx, index=False, na_rep='NA')
    except:
        logger.debug('no excel package installed, skip exporting to excel')
        pass

def prioritize_enhancer(pwout, fn_peak, rep1, rep2, direction, weight, fdr_thres, linkage):
    
    
    # inputs
    # fn_peak  , required if -pri set to 1. ChIP-seq peak file for the regulator of interest in bed format
    # weight, default = 0.5  # Weight to balance the impact of binding and function evidence. The higher the weight, the bigger impact does the binding evidence have (default: 0.5, i.e., binding and functional evidence have equal impact on enhancer prioritization)
    # fdr = 0.05 # Cutoff to select significant transcriptional changes (default: FDR < 0.05). Use Foldchange instead if no FDR is generated (default: Foldchange > 1.5)
    # linkage = 'gb' # gb, pp, pindex  Transcriptional changes to be used for enhancer prioritization, the value can be pp (changes in promoter-proximal region), gb (changes in gene body region), or pindex (changes in pausing index) (default: gb)
    
    # direction = 0  # The expected simultaneous change of expression between enhancers and target genes. 1: the expression changes of enhancers and target genes are in the same direction; -1: the expression changes of enhancers and target genes are in the opposite direction; 0: the expression changes of enhancers and target genes are either in the same or in the opposite direction (default: 0)
    pass

    # variables
    dis_to_p = {}

    if linkage not in {'pp', 'gb', 'pindex'}:
        logger.error(f'invalid linkage specified, valid = pp, gb and pindex')
        return 1


    # format the chromosome of the fn_peak
    fn_peak_new = f'{pwout}/intermediate/peak_file_for_enhancer_prioritization.bed'
    with open(fn_peak) as f, open(fn_peak_new, 'w') as o:
        for i in f:
            chr_, line = i.split('\t', 1)
            chr_ = refine_chr(chr_)
            o.write(f'{chr_}\t{line}')

    # get the nearest peak position for the center
    fn_center = f'{pwout}/eRNA/Enhancer_center.txt'  # previous is intermediate/center.tmp
    fn_enhancer = f'{pwout}/eRNA/Enhancer.txt'
    fn_enhancer_prior = f'{pwout}/eRNA/Enhancer.prioritized.txt'
    
    fn_center_bed = f'{pwout}/intermediate/center.tmp'
    fn_nearest_peak = f'{pwout}/intermediate/nearest.tmp'  # last column is distance

    with open(fn_center) as f, open(fn_center_bed, 'w') as o:
        f.readline()
        for i in f:
            line = i.strip().split('\t')
            chr_, pos, enhancer_id = line[0], line[1], line[3]
            pos = int(pos)
            print(f'{chr_}\t{pos - 1}\t{pos}\t{enhancer_id}', file=o)

    retcode = run_shell(f'bedtools sort -i {fn_center_bed} > {fn_center_bed}.tmp;mv {fn_center_bed}.tmp {fn_center_bed}')
    if retcode:
        logger.error(f'fail to sort file {fn_center_bed}')
    # fn_center format
    # chr	position	FONTOM5	Enhancer_ID
    # 1	737803	N	182
    # 1	2377977	N	126

    cmd = f'bedtools closest -a {fn_center_bed} -b {fn_peak_new} -d > {fn_nearest_peak}'
    logger.debug(cmd)
    retcode = run_shell(cmd)
    with open(fn_nearest_peak) as f:
        for i in f:
            line = i[:-1].split('\t')
            dist = int(line[-1])
            if dist == -1:
                continue
            enhancer_id = line[3]
            if enhancer_id not in dis_to_p:
                dis_to_p[enhancer_id] = dist
            elif dist < dis_to_p[enhancer_id]:
                dis_to_p[enhancer_id] = dist
            
    
    def process_line(line, gn, idx_fc, idx_fdr, prev_max_fc, gmax, tchange, tfdr):
        fc = line[idx_fc]
        gmax_v = None
        if fc not in {'NA', 'Inf', '-Inf'}:
            try:
                fc = float(fc)
                fc_abs = abs(fc)
            except:
                logger.error(f'invalid logFC value found, line = {line}\nexit now...')
                sys.exit(1)
            if fc_abs > prev_max_fc:
                prev_max_fc = fc_abs
        else:
            fc = {'Inf': 1000, '-Inf': -1000, 'NA': 'NA'}[fc]
            gmax_v = {'Inf': 1, '-Inf': -1, 'NA': None}[fc]
            
        if gmax_v is not None:
            gmax[gn] = gmax_v

        if idx_fdr is not None:
            fdr = line[idx_fdr]
            if fdr in {'NA', ''}:
                fdr = 'NA'
            else:
                fdr = float(fdr)
            tfdr[gn] = fdr

        if gn not in tchange:
            tchange[gn] = fc
        elif tchange[gn] == 'NA' or fc_abs > abs(tchange[gn]):
            tchange[gn] = fc
        
        return prev_max_fc
    
    def get_idx(fn, header):
        status = 0
        if 'log2FoldChange' in header:
            idx_fc = header.index('log2FoldChange')
        elif 'log2fc' in header:
            idx_fc = header.index('log2fc')
        else:
            logger.error(f'log2fc column not foun din {fn}, header = \n{header}')
            status = 1
        if rep1 + rep2 == 2:
            idx_fdr = None
        elif 'FDR' in header:
            idx_fdr = header.index('FDR')
        elif 'padj' in header:
            idx_fdr = header.index('padj')
        else:
            logger.error(f'FDR column not found in {fn} header: \n{header}')
            status = 1
        
        return status, idx_fc, idx_fdr

    # no case samples
    if rep2 == 0:
        with open(fn_enhancer) as f, open(fn_enhancer_prior, 'w') as o:
            header = f.readline()[:-1]
            print(f'{header}\tScore', file=o)
            # Enhancer_ID	chr	start	end	center
            for i in f:
                line = i.strip().split('\t', 1)
                enhancer_id = line[0]
                if enhancer_id in dis_to_p:
                    score = 2 / (1 + np.exp(0.0004054651 * dis_to_p[enhancer_id]))
                else:
                    score = 0
                print(f'{i[:-1]}\t{score}', file=o)
        # sort by score, last column
        df = pd.read_csv(fn_enhancer_prior, sep='\t')
        df.sort_values('Score', ascending=False, inplace=True)
        df.to_csv(fn_enhancer_prior, sep='\t', na_rep='NA', index=False)
    else:
        # with case samples
        fn_change = f'{pwout}/known_gene/{linkage}_change.txt'
        fn_enhancer_change = f'{pwout}/eRNA/Enhancer_change.txt'
        # pp, gb, header = 
        # Transcript	Gene	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
        # pindex, header = 
        # Transcript	Gene	log2fc	pvalue	FDR

        prev_max_fc = 0  # prev max abs(log2fc)
        gmax = {}
        tchange = {}
        tfdr = {}
        with open(fn_change) as f:
            header = f.readline()[:-1].split('\t')
            status, idx_fc, idx_fdr = get_idx(fn_change, header)
            if status == 1:
                return 1
            for i in f:
                line = i[:-1].split('\t')
                its, gn = line[:2]
                if not gn:
                    gn = its
                prev_max_fc = process_line(line, gn, idx_fc, idx_fdr, prev_max_fc, gmax, tchange, tfdr)
            for gn, sign in gmax.items():
                gmax[gn] = sign * prev_max_fc
    
    
        # process enhaner_change
        # Enhancer_ID	chr	start	end	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
        temax = 0
        emax = {}
        enhchange = {}
        enhfdr = {}
        with open(fn_enhancer_change) as f:
            header = f.readline()[:-1].split('\t')
            status, idx_fc, idx_fdr = get_idx(fn_enhancer_change, header)
            if status == 1:
                return 1
            for i in f:
                line = i[:-1].split('\t')
                enh_id = line[0]
                temax = process_line(line, enh_id, idx_fc, idx_fdr, temax, emax, enhchange, enhfdr)
        for enh_id, sign in emax.items():
            emax[enh_id] = sign * temax
        
        
        # add bscore and fscore
        n_sam = rep1 + rep2
        if n_sam == 2:
            comp_thres = 0.59
            comp_dict_gene = tchange
            def skip_func(gn):
                if gn == 'NA' or gn not in tchange or comp_dict_gene[gn] == 'NA':
                    return 1
                if abs(comp_dict_gene[gn]) < comp_thres:
                    return 1
                return 0
            def use_default_fscore(gn):
                if enhchange[gn] == 'NA' or (abs(enhchange[gn]) < comp_thres):
                    return 1
                return 0
        else:
            comp_thres = fdr_thres
            comp_dict_gene = tfdr
            def skip_func(gn):
                if gn == 'NA' or gn not in tchange or comp_dict_gene[gn] == 'NA':
                    return 1
                if comp_dict_gene[gn] > comp_thres:
                    return 1
                return 0
            def use_default_fscore(gn):
                if enhfdr[gn] == 'NA' or (enhfdr[gn] > comp_thres):
                    return 1
                return 0

        
        with open(fn_enhancer) as f, open(fn_enhancer_prior, 'w') as o:
            # Enhancer_ID	chr	start	end	center	FANTOM5	associated_gene-FANTOM5	associated_gene-50kb	associated_gene-4DGenome	Cell_Tissue	Detection_method	PMID	closest_gene	distance
            header = f.readline()[:-1]
            print(f'{header}\tFscore\tBscore', file=o)
            for i in f:
                line = i[:-1].split('\t')
                enh_id = line[0]
                bscore = 0
                fscore = 0
                glist_fantom5 = list(set(re.split(r'[,;]', line[6])))
                glist_50k = line[7].split(',')
                glist_4d = line[8].split(',')
                closest_gn = [line[12]]
                
                glist = [[glist_fantom5, 0], [glist_50k, 0], [glist_4d, 0], [closest_gn, 0]]
                
                
                if enh_id in dis_to_p:
                    bscore = 2 / (1 + np.exp(0.0004054651 * dis_to_p[enh_id]))
                    
                if not use_default_fscore(enh_id):
                    for metric in glist:
                        iglist, max_fc = metric
                        for gn in iglist:
                            skip_flag = skip_func(gn)
                            if skip_flag:
                                continue
                            if direction in {1, -1}:
                                concord_enh_gene = enhchange[enh_id] * tchange[gn] * direction
                            else:
                                concord_enh_gene = 1
                            if concord_enh_gene > 0:
                                if abs(tchange[gn]) > abs(metric[1]):
                                    metric[1] = tchange[gn]
                    
                    if direction in {1, -1}:
                        fscore = (enhchange[enh_id] / temax) * sum([_[1] for _ in glist]) / prev_max_fc
                    else:
                        fscore = abs(enhchange[enh_id] / temax) * sum([abs(_[1]) for _ in glist]) / prev_max_fc
                

                print(f'{i[:-1]}\t{fscore}\t{bscore}', file=o)

        df = pd.read_csv(fn_enhancer_prior, sep='\t')
        max_fscore = df['Fscore'].max()
        df['Score'] = df['Bscore'] * weight + df['Fscore']/max_fscore * (1 - weight)
        df.sort_values('Score', ascending=False, inplace=True)
        df.to_csv(fn_enhancer_prior, sep='\t', na_rep='NA', index=False)

    
def show_system_info(home):
    try:
        import getpass
        hostname = str(os.uname())
        # user = os.getlogin()
        user = getpass.getuser()
        # get the group name
        group = os.popen('groups').read().strip()
        # get other basic info
        logger.debug(f'HOME = {home}, hostname = {hostname}, user = {user}, group = {group}')
    except:
        logger.debug('fail to get system info')
    try:
        logger.debug(f'python version = {sys.version}')
        # show python path for executing this script
        logger.debug(f'python path = {sys.executable}')
        logger.info(f'utils path = {__file__}')
        logger.debug(f'current working directory = {os.getcwd()}')
        from __init__ import __version__, __doc__
        logger.info(f'current version = {__version__}, doc = {__doc__}')
    except:
        logger.debug('fail to get python info')
    
def filter_tts_downstream_count(pwout, fn_protein_coding, rep1, rep2, downstream_no_overlap_length, gtf_info_raw):
    downstream_no_overlap_length_str = f'{downstream_no_overlap_length/1000:.0f}k'
    fn_count = f'{pwout}/intermediate/tts_down_{downstream_no_overlap_length_str}.txt'
    fno = f'{pwout}/intermediate/tts_down_{downstream_no_overlap_length_str}_filtered.txt'
    fn_active_gene = f'{pwout}/intermediate/active_gene.txt'
    fno_change = f'{pwout}/known_gene/readthrough_change.txt'
    # header = 
    # Transcript	Gene	chr	strand	last_exon_len	last_exon_s	last_exon_e	last_exon_OmoMYC_1_rRNArm-F4q10	tts_down_50k_OmoMYC_1_rRNArm-F4q10	ratio_OmoMYC_1_rRNArm-F4q10
    # target
    # 1. filter out the protein-coding only
    # 2. for each gene, only keep the active transcript
    # 3. use cmhtest like the pindex change to get the pvalue and odds ratio
    n_info_cols = 7  # first 7 columns are the reference columns
    with open(fn_active_gene) as f:
        active_ts_set = {_.strip().split('\t')[0] for _ in f}
    df = pd.read_csv(fn_count, sep='\t')
    n_orig = len(df)
    df_filter = df.loc[df['Transcript'].isin(active_ts_set)].copy()
    n_filter_active = len(df_filter)
    logger.debug(f'original count = {n_orig}, after filter active = {n_filter_active}, drop = {n_orig - n_filter_active}')
    
    # keep protein coding only
    # proteincoding info get from the ENSEMBL gtf file
    if fn_protein_coding:
        with open(fn_protein_coding) as f:
            f.readline()
            protein_coding_set = {_.split('\t', 1)[0] for _ in f}
        df_filter = df_filter.loc[df_filter['Gene'].isin(protein_coding_set)]
        n_filter_protein_coding = len(df_filter)
        logger.debug(f'after filter protein coding = {n_filter_protein_coding}, drop = {n_filter_active - n_filter_protein_coding}')
    else:
        logger.warning(f'No protein coding gene file provided, skip filtering protein coding genes')

    def get_down_tts_from_gene_info(gene_info):
        start, end, strand, chr_, gn = [gene_info[_] for _ in ['start', 'end', 'strand', 'chr', 'gene_name']]
        s1, e1 = (start, end + downstream_no_overlap_length) if strand == '+' else (start - downstream_no_overlap_length, end)
        return chr_, gn, s1, e1
    ts_without_overlap = get_no_overlap_region_v2(gtf_info_raw, func_get_target_region=get_down_tts_from_gene_info)
    logger.debug(f'transcript without overlap from TSS to {downstream_no_overlap_length_str} downstream of TTS = {len(ts_without_overlap)}')
    df_filter = df_filter.loc[df_filter['Transcript'].isin(ts_without_overlap)]
    n_filter_overlap_genes = len(df_filter)
    logger.warning(f'after filter downstream overlap genes = {n_filter_overlap_genes}, drop = {n_filter_protein_coding - n_filter_overlap_genes}')
    
    df_filter.to_csv(fno, sep='\t', na_rep='NA', index=False)
    
    # run the cmhtest like the pindex change
    cols = list(df_filter.columns)
    idx_last_exon = [i for i, _ in enumerate(cols) if _.startswith('last_exon_') and i >= n_info_cols]
    idx_down = [i for i, _ in enumerate(cols) if _.startswith('tts_down_')]
        
    # df_out = df_filter.iloc[:, :n_info_cols].copy()
    n_sam = rep1 + rep2
    
    if n_sam * 3 + n_info_cols != len(cols):
        logger.error(f'invalid number of columns, expect {n_sam * 3 + n_info_cols}, got {len(cols)}')
        return 1
    
    def cmhtest_for_readthrough(row):
        arr = []
        row = list(row)
        # there should be 3 * n_sam columns
        for i in range(rep1):
            le1, down1 = row[idx_last_exon[i]], row[idx_down[i]]
            for j in range(rep2):
                j += rep1
                try:
                    le2, down2 = row[idx_last_exon[j]], row[idx_down[j]]
                except:
                    logger.info(list(enumerate(cols)))
                    logger.info(f'row cols = {len(row)}, ncol_exp = {len(cols)}')
                    
                    logger.error(f'idx_last_exon = {idx_last_exon}, idx_down = {idx_down}, row = {row}, i = {i}, j = {j}')
                    sys.exit(1)
                if le1 + le2 + down1 + down2 > 0 and down2 * le1 > 0:
                    arr.append([[le2, down2], [le1, down1]])
        odds_ratio, pvalue = cmhtest(arr)
        return np.log2(odds_ratio) if odds_ratio != 'NA' else 'NA', pvalue
    
    
    if n_sam == 2:
        # use fisher exact test
        idx_le1, idx_le2 = idx_last_exon
        idx_down1, idx_down2 = idx_down
        df_filter['pvalue'] = df_filter.apply(lambda x: get_pvalue_2_sample(x, idx_le1, idx_le2, idx_down1, idx_down2), axis=1)
        df_filter = add_FDR_col(df_filter, 'pvalue')
        
        # add odds ratio
        def get_odds_ratio(row):
            le1, le2, down1, down2 = row[idx_le1], row[idx_le2], row[idx_down1], row[idx_down2]
            tmp = le1 * down2
            if all([le1, le2, down1, down2, tmp]):
                return le2 * down1 / tmp
        
        # get the log2fc
        try:
            log2fc = df_filter.apply(lambda x: np.log2(get_odds_ratio(x)), axis=1)
        except:
            e = traceback.format_exc()
            logger.error(e)
            sys.exit(1)
        df_filter['log2fc'] = log2fc
    else:
        # run cmhtest
        df_filter[['log2fc', 'pvalue']] = df_filter.apply(cmhtest_for_readthrough, axis=1, result_type='expand')
        df_filter = add_FDR_col(df_filter, 'pvalue')
        
    df_filter = df_filter.sort_values('FDR')
    
    # reorder the columns, move log3fc, pvalue, FDR before the samples data
    cols_new = cols[:n_info_cols] + ['log2fc', 'pvalue', 'FDR'] + cols[n_info_cols:]
    df_filter = df_filter[cols_new]
    
    df_filter.to_csv(fno_change, sep='\t', na_rep='NA', index=False)
    logger.debug(f'filter readthrough done, output = {fno_change}, nrow = {len(df_filter)}, cols = {list(df_filter.columns)}')

def test_writable(pw):
    try:
        os.makedirs(pw, exist_ok=True)
    except:
        logger.error(f'fail to create {pw}, please check the permission')
        return 1
    fntest = f'{pw}/test_writable.tmp'
    exp_content = 'test1234'
    try:
        with open(fntest, 'w') as f:
            f.write(exp_content)
    except:
        logger.error(f'fail to write to {pw}, please check the permission')
        logger.debug(f'fail to create a file to {fntest}')
        return 1
    with open(fntest) as f:
        readin_content = f.read().strip()
        if readin_content != exp_content:
            logger.error(f'fail to write to {pw}, please check the permission')
            logger.debug(f'exp and actual read in count are different, exp  = {exp_content}, actual = {readin_content}')
            return 1
    os.remove(fntest)
    return 0
def validate_input(args):
    args = vars(args)
    args_filelike = ['in1', 'in2', 'gtf', 'fa', 'pwout', 'pw_bed']
    
    # check if the script is inside of docker
    if not in_docker:
        return 0
    # only need to check absolute path for docker
    err = 0
    for i in args_filelike:
        if i in args and args[i]:
            v = args[i]
            if isinstance(v, str):
                v = [v]
            for j in v:
                if j[0] != '/':
                    logger.error(f'You are running inside of docker, input {i} must be absolute path (starts from root /), got relative path: {j}')
                    err = 1
    return err
