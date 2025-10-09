#! /usr/bin/env python3
import os, sys, hashlib, pickle, gzip, re, json
import sys
def getlogger(logger_name=None, fn_log=None, nocolor=False):
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
        console.setLevel('INFO')
        console.name = 'console'
        logger.addHandler(console)

    if fn_log and 'file' not in handler_names:
        fh_file = logging.FileHandler(fn_log, mode='w', encoding='utf8')
        fh_file.setLevel('DEBUG')
        fh_file.setFormatter(CustomFormatter())
        fh_file.name = 'file'
        logger.addHandler(fh_file)
    return logger
logger = getlogger()

def get_md5_str(s):
    md5 = hashlib.md5()
    md5.update(s.encode())
    return md5.hexdigest()
def refine_chr(chr_):
    chr_ = chr_.strip().lower()
    for _ in ['chromosome', 'chrom', 'chr']:
        chr_ = chr_.replace(_, '')
    return {'23': 'x', '24': 'y', '25': 'mt', 'm': 'mt'}.get(chr_) or chr_

def process_gtf(fn_gtf,  force_rebuild=False, fake_gtf_path=False):
    """
    process the gtf file, get the gene start and end position, and other info
    gtf position is 1-idx, full closed
    fake_gtf_path is just used for building the pre-built pkl file inside of the docker container
    """

    err = {'no_transcript_id': 0, 'no_gene_name': 0, 'invalid_line_format': 0, 'invalid_strand': 0}
    fn_gtf = os.path.realpath(fn_gtf)
    fn_gtf_lb = os.path.basename(fn_gtf).replace('.gz', '').replace('.gtf', '')
    
    if fake_gtf_path is True:
        tmp = os.path.realpath(fn_gtf).split('/')
        fnpure = tmp[-1]
        build = tmp[-2]
        if build.lower() in {'hg19', 'hg38', 'mm10', 'mm39', 'dm3', 'dm6', 'cd10', 'danrer10'}:
            fnpure = f'{build}/{fnpure}'
        gtf_path_in_docker = f'/app/nrsa/ref/{fnpure}'
        fn_gtf_md5 = get_md5_str(gtf_path_in_docker)[:8]
    else:
        fn_gtf_md5 = get_md5_str(fn_gtf)[:8]
    fn_gtf_size = os.path.getsize(fn_gtf)
    fn_gtf_lb = f'{fn_gtf_lb}_md5_{fn_gtf_md5}_size_{fn_gtf_size}'
    logger.info(f'input gtf file = {fn_gtf}, file name string md5sum = {fn_gtf_md5}, file size = {fn_gtf_size}')
    
    logger.info({'fn_gtf': fn_gtf,  'force_rebuild': force_rebuild, 'fn_gtf_lb': fn_gtf_lb})

    gtf_pwout = os.path.dirname(fn_gtf)
    
    if 'gtf' not in fn_gtf.rsplit('.', 2)[-2:]:
        logger.error(f'the gtf file should have the extension of .gtf: {fn_gtf}')
        sys.exit(1)
    fn_gtf_pkl_default = f'{os.path.dirname(fn_gtf)}/{fn_gtf_lb}.gtf_info.pkl'
    fn_gtf_pkl = f'{gtf_pwout}/{fn_gtf_lb}.gtf_info.pkl'
    fn_gtf_meta_json = f'{gtf_pwout}/{fn_gtf_lb}.gtf_meta.json'
    
    if os.path.exists(fn_gtf_pkl_default):
        fn_gtf_pkl = fn_gtf_pkl_default
    logger.info(f'fn_gtf_pkl = {fn_gtf_pkl}, fn_gtf_meta_json = {fn_gtf_meta_json}')
    if os.path.exists(fn_gtf_meta_json):
        # print the content
        with open(fn_gtf_meta_json) as f:
            logger.info(f'gtf_meta_json = \n{f.read()}')

    if os.path.exists(fn_gtf_pkl) and not force_rebuild:
        logger.info(f'loading gtf from pickle: {fn_gtf_pkl}')
        with open(fn_gtf_pkl, 'rb') as f:
            gtf_info = pickle.load(f)
        return gtf_info, err
    elif os.path.exists(fn_gtf_pkl) and force_rebuild:
        logger.warning(f'force rebuild the gtf pkl file, ignore existing pkl file: {fn_gtf_pkl}')
    
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
                    logger.info(f'error during processing gtf {line_err}: {i}')
                
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
    logger.info(f'g@merged {n_merged} transcripts with same start and end position')
    # logger.info(meta)
    
    with open(fn_gtf_pkl, 'wb') as o:
        logger.info(f'dump to {fn_gtf_pkl}')
        pickle.dump(res, o)
    
    with open(fn_gtf_meta_json, 'w') as o:
        logger.info(f'dump to {fn_gtf_meta_json}')
        json.dump(meta, o, indent=4)


    err_total = sum(err.values())
    if err_total:
        logger.info(f'error in parsing gtf file: {err}')
    
    return res, err



if __name__ == "__main__":
    import argparse as arg
    from argparse import RawTextHelpFormatter
    ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    ps.add_argument('fn', help="""file name for gtf file""")
    ps.add_argument('-force', '-f', help="""force rebiuld the pkl file""", action='store_true')
    
    args = ps.parse_args()

    process_gtf(args.fn, force_rebuild=args.force, fake_gtf_path=True)
