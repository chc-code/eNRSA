#! /usr/bin/env python3
"""
process PROseq data
"""
def getargs():
    import argparse as arg
    from argparse import RawTextHelpFormatter
    ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    ps.add_argument('-in1', help="""required if -design_table is not specified, read alignment files in bed (6 columns) or bam format for condition1, separated by space""", nargs='+')
    ps.add_argument('-in2', help="""read alignment files in bed (6 columns) or bam format for condition2, separated by space""", nargs='*')
    ps.add_argument('-pwout', '-output', help="""required, output/work directory""", required=True)
    ps.add_argument('-organism', '-org',  help="""required. define the genome: built-in = hg19, hg38, mm10, mm39, dm3, dm6, ce10, or danRer10. if other organism used, must also specify the -gtf and -fa""", required=True)

    ps.add_argument('-design_table', '-design',  help="""Optional, desgin table in tsv format. 2 sections, first section is the sample information, 2 or 3 columns. col1=bed/bam full file path, col2 = group name, col3=optional, the batch_group. Second section is the comparison information, 2 columns, col1 should start with @@ used to identify this line as comparison definition. col1=group name used as case, col2 = group name used as control. e.g. @@case\\tcontrol. If only need to process a single group, use @@null\\tgroup_name""", nargs='?')
    ps.add_argument('-batches', '-batch', '-b',  help="""Optional, if -in1 / -in2 is specified, and the need to adjust the batch effect, add the batch label for each sample. the total argument number should match with in1 + in2, sep = space""", nargs='*')

    ps.add_argument('-gtf', help="""user specified GTF file, if not specified, will use the default GTF file for the organism""")
    ps.add_argument('-fa', help="""Full path for the fasta file for the genome. If not specified, will search under the fa folder under the package directory. e.g. organism is hg19, then the default fasta file should be fa/hg19.fa""")
    ps.add_argument('-f1', help="""normalization factors for samples of condition1, separated by space. If no normalization factor is specified, we will use DEseq2 default method to normalize""", nargs='*', type=float)
    ps.add_argument('-f2', help="""normalization factors for samples of condition2, same as -f1""", nargs='*', type=float)
    ps.add_argument('-up', '-pp_up', help="""define the upstream of TSS as promoter (bp, default: 500)""", type=int, default=500)
    ps.add_argument('-down', '-pp_down', help="""define the downstream of TSS as promoter (bp, default: 500)""", type=int, default=500)
    ps.add_argument('-gb', '-gb_start', help="""define the start of gene body density calculation (bp, default: 1000)""", type=int, default=1000)
    ps.add_argument('-tts_padding', help="""define the downstream distance of TTS to get the TTS region count (bp, default: 2000)""", default=2000, type=int)
    ps.add_argument('-min_gene_len', help="""define the minimum length of a gene to perform the analysis (bp, default: 1000). Genes with length less than this will be ignored""", type=int, default=1000)
    ps.add_argument('-window', '-window_size', help="""define the window size (bp, default: 50)""", type=int, default=50)
    ps.add_argument('-step', '-step_size', help="""define the step size (bp, default: 5)""", type=int, default=5)
    ps.add_argument('-bench', '-test', help="""bench mode, only process the first 10 genes for testing purpose""", action='store_true')
    ps.add_argument('-force', help="""ignore the previous pre-count result for the bed files""", action='store_true')
    # ps.add_argument('-overwrite', help="""if the mapped reads count file already exists, overwrite it, default is reuse the old ones""", action='store_true')
    ps.add_argument('-verbose', '-v', help="""verbose mode, print more information""", action='store_true')
    ps.add_argument('-demo', help="""demo mode, skip running the get_mapped_reads step if already exist""", action='store_true')
    ps.add_argument('-ignore', help="""ignore the existing pre-counting results""", action='store_true')
    ps.add_argument('-sorted', help="""the input bed files are sorted, skip the sorting step""", action='store_true')
    # ps.add_argument('-testfunc', help="""test the new functions,debug mode""", action='store_true')
    ps.add_argument('-tts_down', help="""TTS downstream length for detecting readthrough""", type=int, default=50000)
    args = ps.parse_args()
    
    return args

args = getargs()  # just to prevent loading below utils if just calling help

import time
s = time.time()
import sys, os, re
import logging
import pickle
import json
import traceback
import gc
sys.dont_write_bytecode = True

from utils import check_dependency, build_idx_for_fa,  process_gtf,  change_pp_gb, change_pindex, draw_box_plot, draw_heatmap_pindex, draw_heatmap_pp_change, get_FDR_per_sample, add_value_to_gtf, time_cost_util, parse_design_table, get_alternative_isoform_across_conditions, build_design_table, show_system_info, filter_tts_downstream_count, test_writable, validate_input

from utils import Analysis, process_bed_files


time_cost = {}
now = time.time()

s = now

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

def updatelogger(logger, fn_log, terminal_level=None):
    handlers = list(logger.handlers)
    handler_name_d = {v.name: idx for idx, v in enumerate(handlers)}
    fh_console = handlers[handler_name_d['console']]
    formatter = fh_console.formatter
    valid_levels = {'DEBUG', 'INFO', 'WARNING', 'ERROR', 'FATAL'}
    if terminal_level:
        terminal_level = terminal_level.upper()
        if terminal_level not in valid_levels:
            logger.error(f'invalid logging level: {terminal_level} ')
            return logger
    if fn_log and 'file' not in handler_name_d:
        fh_file = logging.FileHandler(fn_log, mode='w', encoding='utf8')
        fh_file.setLevel('DEBUG')
        fh_file.setFormatter(formatter)
        fh_file.name = 'file'
        logger.addHandler(fh_file)
    if terminal_level is not None:
        fh_console.setLevel(terminal_level)
    return logger

logger = getlogger(logger_name='NRSA')


def bench(s, lb):
    time_cost.setdefault(lb, 0)
    now = time.time()
    time_cost[lb] += now - s
    return now


def collect_previous_count(analysis, fls):
    """
    get the previous done _raw, _FDR file
    """
    pp_str = {} # key = transcript_id, v = the peak count in promoter region for each file. 
    gb_str = {} # key = transcript_id, v = [the peak count in gene body region, ratio] for each file.

    for fn_lb, _ in fls:
        fh_bed_peaks = analysis.out_fls['bed_peaks'][fn_lb]['fh']
        # Transcript	Gene	ppc	ppm	ppd	pps	gbc	gbm	gbd	pauseIndex
        fh_bed_peaks.seek(0)
        header = fh_bed_peaks.readline()[:-1].split('\t') # header
        idx_ppc = header.index('ppc')
        for line in fh_bed_peaks:
            row = line.strip().split('\t')
            transcript_id = row[0]
            ppc, gbc, gbd = [row[_ + idx_ppc] for _ in [0, 4, 6]]
            pp_str.setdefault(transcript_id, [])
            pp_str[transcript_id].append(ppc)
            gb_str.setdefault(transcript_id, [])
            gb_str[transcript_id] += [gbc, gbd]

    return pp_str, gb_str


def main(args):
    """
    args: an object with attributes equiv to the argparse object, must specify 
    """

    # check if the attributes are present
    defined_attrs = vars(args)
    required_attrs = {'in1', 'pwout', 'organism', 'pw_bed'}
    
    optional_attrs = {
        'in2': None,
        'gtf': None,
        'fa': None,
        'f1': None,
        'f2': None,
        'up': 500,
        'down': 500,
        'gb': 1000,
        'min_gene_len': 1000,
        'window': 50,
        'step': 5,
        'bench': False,
        'batches': None,
        'ignore': False,
        'tts_padding': 2000,
    }
    pw_bed = args.pw_bed
    pwout = args.pwout
    
    status = validate_input(args)
    if status:
        logger.error(f"Invalid input arguments")
        sys.exit(1)
    
    # test writtable
    write_error = test_writable(pwout)
    if write_error:
        logger.error(f"Output directory is not writable: {pwout}")
        sys.exit(1)
    
    
    missing = required_attrs - set(defined_attrs)
    exist = required_attrs & set(defined_attrs)
    for attr in exist:
        if getattr(args, attr) is None:
            missing.add(attr)
    if missing:
        logger.error(f"Missing required attributes: {sorted(missing)}")
        return
    for attr, default in optional_attrs.items():
        if attr not in defined_attrs:
            setattr(args, attr, default)
    
    logger.debug(f'g@args={vars(args)}')
    
    benchmode = args.bench
    # check dependencies
    dependency_status = check_dependency()
    if dependency_status:
        logger.error("Dependency check failed")
        sys.exit(1)
    demo = vars(args).get('demo', False)

    for lb, d in zip(['Output', 'Intermediate', 'Known Genes'], [pwout, f'{pwout}/intermediate', f'{pwout}/known_gene']):
        if not os.path.exists(d):
            logger.debug(f'Making {lb} directory: {d}')
            os.makedirs(d)


    build_design_table(pwout, args.in1, args.in2, args.batches)

    analysis = Analysis(args, skip_get_mapped_reads=demo)
    if analysis.status:
        logger.error("Exit now")
        sys.exit(1)
    rep1 = len(analysis.control_bed)
    rep2 = len(analysis.case_bed) if analysis.case_bed else 0
    window_size = analysis.config['window_size']
    factor1 = analysis.config['normalize_factor1']
    factor2 = analysis.config['normalize_factor2']
    factor_flag = 0 if factor1 is None else 1

    pro_up, pro_down, gb_down_distance, min_gene_len, tts_padding = [analysis.config[_] for _ in ['pro_up', 'pro_down', 'gb_start', 'min_gene_len', 'tts_padding']]

    # prepare fasta file
    fn_fa = analysis.ref['fa']
    fa_idx = build_idx_for_fa(fn_fa)
    fh_fa = open(fn_fa, 'r')


    # process the GTF file
    logger.info(f"Processing GTF file: {analysis.ref['gtf']}")
    gtf_info, fn_tss, fn_tss_tts, err = process_gtf(analysis.ref['gtf'], pwout=analysis.pwout)
    gtf_info_raw = gtf_info
    err_total = sum(err.values())
    if err_total:
        logger.warning(f'Error encountered while processing GTF file')
        print(json.dumps(err, indent=4))
    len1 = len(gtf_info)
    if len1 == 0:
        logger.error("No gene information found in the GTF file")
        sys.exit(1)
        
    # filter out the genes with length less than the minimum gene length
    gtf_info_new = {}
    short_genes = 0
    merged_transcripts = 0

    tmp = sorted(gtf_info)
    skipped_ts = {'invalid_chr': 0,  'short': 0}
    # ts_prefix = {}
    for k in tmp:
        v = gtf_info[k]
        if fa_idx and v['chr'] not in fa_idx:
            skipped_ts['invalid_chr'] += 1
            continue
        if v['end'] - v['start'] + 1 > min_gene_len: # not include if len is 1000 (min_gene_len)
            merged_transcripts += k.count(';')
            gtf_info_new[k] = add_value_to_gtf(v, pro_up, pro_down, gb_down_distance, tts_padding) 
        else:
            skipped_ts['short'] += 1
            short_genes += 1
    total_skipped = sum(skipped_ts.values())
    len1 = len(gtf_info_new)
    logger.debug(f'initial gtf dict = {len(gtf_info)}, total skipped = {total_skipped}, detail = {skipped_ts} merged transcript = {merged_transcripts}, total count before merge = {len1 + merged_transcripts}, current ts count = {len1}')
    gtf_info = gtf_info_new
    
    fls = analysis.input_fls # element is [fn_lb, fn_bed]
    sam_order = [_[0] for _ in fls]
    n_gene_cols = analysis.n_gene_cols
    
    reuse_pre_count = not vars(args).get('ignore', False)

    fn_count_pp_gb = analysis.out_fls['count_pp_gb']
    line_count = 0
    if os.path.exists(fn_count_pp_gb):
        with open(fn_count_pp_gb) as f:
            for _ in f:
                line_count += 1

    # logger.warning('modify here')
    # fn_protein_coding = analysis.ref['protein_coding']
    # filter_tts_downstream_count(pwout, fn_protein_coding, rep1, rep2)
    # sys.exit(1)
    
    if not (line_count > 10 and demo):
        logger.info(f'Getting pp_gb count')
        tts_down = args.tts_down
        logger.info(f'TTS downstream length for detecting readthrough = {tts_down/1000} k')
        pp_str, gb_str = process_bed_files(analysis, fls, gtf_info, gtf_info_raw, fa_idx, fh_fa, reuse_pre_count=reuse_pre_count, downstream_no_overlap_length=tts_down)

        # close file handle
        for fn_lb, fn_bed in fls:
            analysis.out_fls['bed_peaks'][fn_lb]['fh'].close()
        fh_fa.close()
        # calculate FDR
        logger.info('Calculating FDR')
        pause_index_str = {}
        for fn_lb, fn_bed in fls:
            fn_peak = analysis.out_fls['bed_peaks'][fn_lb]['fn']
            fn_fdr = analysis.out_fls['fdr'][fn_lb]['fn']
            pause_index_str = get_FDR_per_sample(fn_lb, fn_peak, fn_fdr,  pause_index_str)
    else:
        logger.debug('Re-use previous mapped reads count')
        pp_str, gb_str = collect_previous_count(analysis, fls)
    
        logger.debug('collecting pause index from FDR files')
        pause_index_str = {}  # key = transcript_id
        for fn_lb, _ in fls:
            fn_peak = analysis.out_fls['bed_peaks'][fn_lb]['fn']
            fn_fdr = analysis.out_fls['fdr'][fn_lb]['fn']
            if not os.path.exists(fn_fdr):
                pause_index_str = get_FDR_per_sample(fn_lb, fn_peak, fn_fdr,  pause_index_str)
            with open(fn_fdr) as f:
                header = f.readline().strip().split('\t')
                idx_pindex, idx_pval, idx_fdr = header.index('pauseIndex'), header.index('pvalue'), header.index('FDR')
                for i in f:
                    line = i[:-1].split('\t')
                    transcript_id, pindex, pval, fdr = [line[_] for _ in [0, idx_pindex, idx_pval, idx_fdr]]
                    pause_index_str.setdefault(transcript_id, {})[fn_lb] = [pindex, pval, fdr]
        

    pause_index_str_new = {}
    for transcript_id, v1 in pause_index_str.items():
        v2 = []
        for isam in sam_order:
            v3 = v1.get(isam, ['NA', 'NA', 'NA'])
            v2 += v3
        pause_index_str_new[transcript_id] = v2
    del pause_index_str
    gc.collect()
    pause_index_str = pause_index_str_new
    
 
    # count_pp_gb
    logger.info('Building count_pp_gb')
    prefix = 'longeRNA-' if analysis.longerna else ''
    header = analysis.gene_cols + ['chr', 'start', 'end', 'strand']
    header_pp = []
    header_gb = []
    for fn_lb, fn_bed in fls:
        header_gb += [f'gbc_{fn_lb}', f'gbd_{fn_lb}']
        header_pp.append(f'ppc_{fn_lb}')
    header += header_pp + header_gb

    
    with open(fn_count_pp_gb, 'w') as o:
        print('\t'.join(header), file=o)
        for transcript_id in pp_str:
            gene_info = gtf_info[transcript_id]
            # row = []
            row = [transcript_id]
            if not analysis.longerna:
                row.append(gene_info['gene_name'])
            row +=  [gene_info['chr'], str(gene_info['start']), str(gene_info['end']), gene_info['strand']]
            row += pp_str[transcript_id] + gb_str[transcript_id]
            print('\t'.join(row), file=o)


    fn_count_pp_gb_debug = f'{pwout}/intermediate/count_pp_gb_debug.txt'
    with open(fn_count_pp_gb_debug, 'w') as o:
        header = analysis.gene_cols + ['chr', 'start', 'end', 'strand', 'pp_range', 'gb_range'] + header_pp + header_gb
        print('\t'.join(header), file=o)
        for transcript_id in pp_str:
            gene_info = gtf_info[transcript_id]
            pp_range = f'{gene_info["chr"]}:{gene_info["pp_start"]}-{gene_info["pp_end"]}:{gene_info["strand"]}'
            gb_range = f'{gene_info["chr"]}:{gene_info["gb_start"]}-{gene_info["gb_end"]}:{gene_info["strand"]}'
            row = [transcript_id, gene_info['gene_name'], gene_info['chr'], str(gene_info['start']), str(gene_info['end']), gene_info['strand'], pp_range, gb_range]
            row += pp_str[transcript_id] + gb_str[transcript_id]
            print('\t'.join(row), file=o)

    del pp_str
    del gb_str
    gc.collect()

    # if time_cost:
    #     tmp = json.dumps(time_cost, indent=3)
    #     print(tmp)
    
    fno_prefix = 'longeRNA-' if analysis.longerna else ''

    fno_ppchange = f'{pwout}/known_gene/pp_change.txt'
    fno_gbchange = f'{pwout}/known_gene/gb_change.txt'
    fn_gbc_sum = f'{pwout}/intermediate/gbc_sum.json'
    if not (os.path.exists(fno_ppchange) and os.path.exists(fno_gbchange) and os.path.exists(fn_gbc_sum) and demo):
        logger.info('Change_pp_gb')
        change_pp_gb(n_gene_cols, fn_count_pp_gb, analysis.pwout, rep1, rep2, window_size, factor1=factor1, factor2=factor2, factor_flag=factor_flag, islongerna=analysis.longerna)
    else:
        logger.debug(f'skip change_pp_gb due to demo mode')
    
    logger.debug('Dump pindex.txt')
    header_extra = []
    for fn_lb, fn_bed in fls:
        header_extra += [f'{fn_lb}-pindex', f'{fn_lb}-pvalue', f'{fn_lb}-FDR']
    header = analysis.gene_cols + header_extra
    fno =  os.path.join(analysis.known_gene_dir, fno_prefix + 'pindex.txt') 
    fno_pindex_change =  os.path.join(analysis.known_gene_dir, fno_prefix + 'pindex_change.txt') 
    with open(fno, 'w') as o:
        print('\t'.join(header), file=o)
        for transcript_id, data in pause_index_str.items():
            row = [transcript_id]
            if not analysis.longerna:
                row.append(gtf_info[transcript_id]['gene_name'])
            row += data
            print('\t'.join(map(str, row)), file=o)

    # boxplot
    logger.info(f'plotting boxplot')
    draw_box_plot(n_gene_cols, analysis.pwout, 'boxplot', rep1, rep2)
    if rep2 > 0:
        # change_pindex
        logger.debug('Running change_pindex')
        change_pindex(fno_prefix, n_gene_cols, fn_count_pp_gb, fno_pindex_change, rep1, rep2, window_size, factor1=None, factor2=None, factor_flag=0)
        
        
        # simple_heatmap
        logger.info(f'plotting pindex heatmap')
        draw_heatmap_pindex(analysis.pwout)
        

        # heatmap
        logger.info(f'plotting heatmap for pp_change')
        # "perl heatmap.pl -w $pwout -i $list -in1 $cond1_str -in2 $cond2_str -m $genome -n $tname";
        fno = f'{pwout}/known_gene/heatmap.pdf'
        if not (os.path.exists(fno) and demo):
            draw_heatmap_pp_change(n_gene_cols, analysis.pwout, pw_bed,  fls_ctrl=analysis.control_bed, fls_case=analysis.case_bed, fn_tss=fn_tss, region_size=5000, bin_size=200, outname='heatmap', skipe_bedtools_coverage=demo)
        else:
            logger.debug(f'skip plot heatmap due to demo mode')

    tmp = json.dumps(time_cost_util, indent=4)
    logger.info(tmp)

    # # get alternative isoforms
    logger.info('Getting alternative TSS isoforms')
     
    fn_count_pp_gb = f'{pwout}/intermediate/count_pp_gb.txt'
    fn_count_tts = f'{pwout}/intermediate/count_tts.filtered.txt'
    # logger.debug(f'pwout = {pwout}')
    # logger.debug(vars(analysis))
    get_alternative_isoform_across_conditions(fn_count_pp_gb, pwout, pw_bed, rep1, rep2)
    get_alternative_isoform_across_conditions(fn_count_tts, pwout, pw_bed, rep1, rep2, tts_padding=tts_padding)

        # filter and add FDR
    if rep2 > 0:
        logger.debug(f'filtering on TTS downstream readthrough table')
        downstream_no_overlap_length = args.tts_down
        fn_protein_coding = analysis.ref['protein_coding']
        filter_tts_downstream_count(pwout, fn_protein_coding, rep1, rep2, downstream_no_overlap_length, gtf_info_raw)
    

    
if __name__ == "__main__":
    args = getargs()
    pwout = args.pwout
    os.makedirs(pwout, exist_ok=True)

    fn_log = f'{pwout}/NRSA.run.log'
    fn_log_base = os.path.basename(fn_log)
    terminal_level = 'DEBUG' if args.verbose else None
    logger = updatelogger(logger, fn_log, terminal_level=terminal_level)

    if os.path.exists(fn_log_base):
        os.unlink(fn_log_base)
    if os.path.exists(fn_log):
        os.symlink(fn_log, fn_log_base)

    show_system_info()
    logger.debug(f'working in {os.getcwd()}')
    logger.debug(f'inpu args = {vars(args)}')
    
    pwout_raw = os.path.realpath(args.pwout)
    logger.debug(f'pw_out_raw = {pwout_raw}')
    args.pwout = pwout_raw
    args.pwout_raw = pwout_raw
    args.pw_bed = f'{pwout_raw}/bed'
    if not os.path.exists(args.pw_bed):
        os.makedirs(args.pw_bed, exist_ok=True)

    arg_list = parse_design_table(args)
    retcode = 0
    for comp_str, iargs in arg_list:
        if comp_str:
            logger.info(f'g@now running {comp_str}')
        try:
            retcode = main(iargs) or retcode
        except:
            e = traceback.format_exc()
            if 'sys.exit(1)' not in e:
                logger.error(f'Error found during running: log file= {fn_log}, error info = \n\n{e}')
            else:
                logger.error(f'exit now ...')
            sys.exit(1)

    if not retcode:
        logger.debug(f'g@script finished without error')
    else:
        logger.debug(f'error encountered during running')

    sys.exit(retcode)
