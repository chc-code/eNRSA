#! /usr/bin/env python3
"""
get the reads count of specified region(s) from the bed/bam file
"""
import argparse as arg
from argparse import RawTextHelpFormatter
ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
ps.add_argument('fls', help="""read alignment file(s) in bam/bed format, can be multiple, sep by space, if bam is specified, will be onverted to bed first""", nargs='+')
ps.add_argument('--region', '-r', help="""region(s) to count, in the format of chr:start-end:strand, can be multiple, sep by space, coordinate is 1-indexed""", nargs='+')
ps.add_argument('-fn_region', '-R', '-bed', help="""bed file containing the region(s) to count, in the format of chr start end, strand. optional, 5th column can be the region name. 0-indexed, will be ignored if -r is specified""")
ps.add_argument('--verbose', '-v', help="""show debug level log""", action='store_true')
ps.add_argument('-demo', '-validate', '-test', help="""get the count of the region manually, used to validate the result, only accept -region as input""", action='store_true')
args = ps.parse_args()


import os, sys, re
import gzip
import logging
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))

from nrsa.utils import get_peak_method1, get_peak_method2, pre_count_for_bed, refine_chr, process_input


def getlogger(fn_log=None, logger_name=None, nocolor=False, debug=False):
    
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
        console_level = 'DEBUG' if debug else 'INFO'
        console.setLevel(console_level)
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
logger = getlogger(fn_log='get_region_count.log', logger_name='NRSA')
if args.verbose:
    logger = updatelogger(logger, None, 'DEBUG')
logger.debug(vars(args))

def get_region(regions_terminal, regions_bed):
    regions = {}
    invalid_input = []
    
    def validate(orig_str, chr_, start, end, strand):
        chr_ = refine_chr(chr_)
        if strand not in {'+', '-'}:
            invalid_input.append(f'invalid_strand\t{orig_str}')
            return 1
        try:
            start, end = int(start), int(end)
        except:
            invalid_input.append(f'start/end_not_int\t{orig_str}')
            return 1
        
        if start > end:
            invalid_input.append(f'start_larger_than_end\t{orig_str}')
            return 1
        
        return chr_, start, end, strand

    
    if regions_terminal:
        for i in regions_terminal:
            i = re.sub(r'[\W_]+', ':', i[:-1]) + i[-1]
            tmp = re.split(r':', i)
            if len(tmp) != 4:
                invalid_input.append(f'missing_fields\t{i}')
                continue
            tmp1 = validate(i, *tmp)
            if tmp1 == 1:
                continue
            chr_, start, end, strand = tmp1
            k = f'{chr_}:{start}-{end}:{strand}'
            regions[k] = [chr_, start, end, strand]
    elif regions_bed:
        with open(regions_bed) as f:
            for line in f:
                line = line.strip()
                tmp = re.split(r'\s+', line)
                if len(tmp) < 4:
                    invalid_input.append(f'missing_fields\t{line}')
                    continue
                tmp1 = validate(line, *tmp[:4])
                if tmp1 == 1:
                    continue
                chr_, start, end, strand = tmp1
                start += 1
                k = tmp[4] if len(tmp) > 4 else f'{chr_}:{start}-{end}:{strand}'
                regions[k] = [chr_, start, end, strand]
    if len(invalid_input) > 0:
        logger.error(f'Invalid input:\n\t{"\n\t".join(invalid_input)}')

    return regions

def main(pwd, fls, regions):
    
    # get the pre-count for each bed file
    pw_bed = f'{pwd}/bed'
    if not os.path.exists(pw_bed):
        os.makedirs(pw_bed, exist_ok=True)
    
    res = {k: [0 for _ in range(len(fls))] for k in regions} # k = region_id, v = count in each file
    for i_file, (fn_lb, fn_out_bed) in enumerate(fls):
        count_per_base, count_bin = pre_count_for_bed(fn_lb, fn_out_bed, pw_bed, bin_size=200, reuse=True)
        bin_size = count_bin['bin_size']
        for k, [chr_, start, end, strand] in regions.items():
            # get_peak_method1(count_per_base, chr_, strand_idx, s, e)
            # get_peak_method2(count_per_base, count_bin, chr_, strand_idx, s, e, bin_size)
            if end - start > bin_size:
                strand_idx = 0 if strand == '+' else 1
                res[k][i_file] = get_peak_method2(count_per_base, count_bin, chr_, strand_idx, start, end, bin_size)
            else:
                res[k][i_file] = get_peak_method1(count_per_base, chr_, strand_idx, start, end)


    # dump the results to tsv, row = region col = fn_lb
    fn_res = f'{pwd}/region_count.tsv'
    fls_lb = [_[0] for _ in fls]
    header = f'region\t{"\t".join(fls_lb)}'
    with open(fn_res, 'w') as f:
        print(header, file=f)
        for k, v in res.items():
            f.write(f'{k}\t{"\t".join(map(str, v))}\n')
    
    # if the result is very small, < 10, print to termnial
    if len(res) < 10:
        # get the max length of each column
        col_width = [0 for _ in range(len(fls_lb) + 1)]
        header_l = ['region'] + fls_lb
        for i, k in enumerate(header_l):
            if len(k) > col_width[i]:
                col_width[i] = len(k)
        
        for k, v in res.items():
            for i, j in enumerate([k] + v):
                if len(str(j)) > col_width[i]:
                    col_width[i] = len(str(j))
    
        # print the result
        extra_space = 3
        def format_line(l, col_width):
            return ' '.join([f'{i:<{j+extra_space}}' for i, j in zip(l, col_width)])
        print(format_line(header_l, col_width))
        for k, v in res.items():
            print(format_line([k] + v, col_width))


def get_count_manual(fn_bed, chr_, start, end, strand):
    """
    get the count by iterate over the bed file, used to validat the result
    """
    # example:
        # 1	235824344	236029220	-	LYST	(gb_start, gb_end)
        # count of 4 files (hg19 , 4104-AW-1 to 4): 1934	1686	1986	1327
    
    with gzip.open(fn_bed, 'rt') if fn_bed.endswith('.gz') else open(fn_bed) as f:
        ct = 0
        chr_found = 0
        
        for line in f:
            line = line.strip().split('\t')
            chr_ = refine_chr(line[0])
            if chr_ != chr_:
                if chr_found:
                    break
                continue
            if strand != line[5]:
                continue
            chr_found = 1
            s, e = int(line[1]), int(line[2])
            read_end = e if strand == '+' else s
            if start <= read_end <= end:
                ct += 1
            if s > end:
                break
        
    return ct
            


if __name__ == "__main__":
    pwd = os.getcwd()
    fls_raw = args.fls
    regions_terminal = args.region
    regions_bed = args.fn_region
    
    if args.demo and not regions_terminal:
        logger.error(f'--demo can only be used with --region')
        sys.exit(1)
    if args.demo:
        region_bed = None
    
    
    logger.info('Parsing regino input')
    regions = get_region(regions_terminal, regions_bed)
    if len(regions) == 0:
        logger.error(f'No region specified')
        sys.exit(1)
    
    fls = []
    fls_dedup = set()
    for fn in fls_raw:
        fn = os.path.realpath(fn)
        if fn in fls_dedup:
            continue
        fls_dedup.add(fn)
        fls.append(fn)
    # logger.debug(fls)
    # sys.exit(1)
    
    # convert bam to bed if needed
    logger.info('Processing input files')
    fls = process_input(pwd, fls, respect_sorted=True)  # each element is [fn_lb, fn_bed]
    if fls is None:
        logger.error(f'No valid input file')
        sys.exit(1)
    if not args.demo:
        main(pwd, fls, regions)
    else:
        logger.warning(f'demo mode, will get the count manually')
        res = {}
        for fn_lb, fn_bed in fls:
            for region_name, [chr_, start, end, strand] in regions.items():
                ict = get_count_manual(fn_bed, chr_, start, end, strand)
                res.setdefault(region_name, []).append(ict)
        
        print([_[0] for _ in fls])
        for k, v in res.items():
            print(k, v)
