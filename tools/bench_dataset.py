#! /usr/bin/env python3
"""
benchmark for different datasets, run in cqs3
"""
import sys
import os, re
import json
import traceback # use trackback.format_exc() to capture info
home = os.path.expanduser("~")


hostname = os.uname()[1].split('.')[0]

import sys
def getlogger(fn_log=None, logger_name=None, nocolor=False):
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

def prepare_config(pwd):
    pw_gtf = f'{pwd}/gtf'
    pw_bed = f'{pwd}/bed'
    pw_sh = f'{pwd}/sh'
    pw_out = f'{pwd}/out'
    flist_other_pw = '/nobackup/tansey_lab/from-nobackup-h_vangard_1-wangj52-Bill/Andrea-G401-KYM-PROseq/bowtie2-recover/result/'
    
    for ipw in [pw_gtf, pw_bed, pw_sh, pw_out]:
        os.makedirs(ipw, exist_ok=True)

    flist_size = [
        '/workspace/wangj52/PROseq-data/AML-GSE83660/GSM2212051_sorted_rRNArm.bam',
        '/workspace/wangj52/PROseq-data/Bcell-GM12004-GSM980646/bams/GSM980646-GM12004_PROseq_sorted_rRNArm-F4q10.bam',
        '/workspace/wangj52/PROseq-data/CD4-GSM1613181/bams/CD4-PRO-UT.rRNArm.F4q10.sorted.bam',
        '/workspace/wangj52/PROseq-data/Daudi-Pankaj/bams/Daudi-PRO-DMSO.rRNArm.F4q10.sorted.bam',
        '/workspace/wangj52/PROseq-data/K562-GSM1480327/bams/K562-PRO-UT.rRNArm.F4q10.sorted.bam',
    ]
    flist_size = [['size', fn] for fn in flist_size]
    
    flist_other = [
        ['ctrl', '4104-AW-1_sorted.bam'],
        ['ctrl', '4104-AW-3_sorted.bam'],
        ['case', '4104-AW-2_sorted.bam'],
        ['case', '4104-AW-4_sorted.bam'],
    ]
    flist_other = [[lb, f'{flist_other_pw}/{fn}'] for lb, fn in flist_other]
    flist = flist_size + flist_other
    
    # the datasets for testing file size
    file_size_json = f'{pwd}/bam_file_size.json'
    if os.path.exists(file_size_json):
        with open(file_size_json, 'r') as f:
            file_size = json.load(f)
    else:
        file_size = {}
        
    file_not_exist = []
    files_for_size = []
    files_case_ctrl = {}

    force_rerun_get_size = False
    for lb, fn in flist:
        fn_bed = f'{pw_bed}/{os.path.basename(fn).replace(".bam", ".bed")}'
        fn_bed = re.sub(r'[_.-]?sort(?:ed)?', '', fn_bed)
        fn_bed = re.sub(r'[_.\-](rRNArm|F4q10)', '', fn_bed)

        if not os.path.exists(fn) and not os.path.exists(fn_bed):
            file_not_exist.append(fn)
            continue


        if not os.path.exists(fn_bed):
            logger.warning(f'converting bam to bed: {lb} - {os.path.basename(fn)}')
            os.system(f'bedtools bamtobed -i {fn} |bedtools sort > {fn_bed}')   


        if fn_bed not in file_size or force_rerun_get_size:
            if not os.path.exists(fn_bed):
                file_not_exist.append(fn_bed)
                continue
            isize = os.path.getsize(fn_bed)
            # convert to GB if > 1GB else convert to MB, round to 2 decimal, also include the unit
            isize = f'{isize/1024/1024/1024:.1f} GB' if isize > 1024*1024*1024 else f'{isize/1024/1024:.1f} MB'
            # get the reads count using samtools
            # cmd = f'samtools view -c {fn}'
            cmd = f'cat {fn_bed}|wc -l'
            logger.info(f'getting reads count: {lb} - {os.path.basename(fn_bed)}')
            read_count = int(os.popen(cmd).read().strip()) / 1000/1000 # convert to million
            file_size[fn_bed] = {'size': isize, 'reads': f'{read_count:.1f}M'}
        else:
            isize = file_size[fn_bed]['size']
        
        if lb == 'size':
            files_for_size.append([isize, fn_bed])
        else:
            files_case_ctrl.setdefault(lb, []).append(fn_bed)

    # logger.info(files_case_ctrl)
    if file_not_exist:
        logger.error(f'files not exist n = {len(file_not_exist)} : {file_not_exist}')
        sys.exit(1)
    with open(file_size_json, 'w') as f:
        json.dump(file_size, f, indent=4)
    
    # build the -in1 and -in2, -gtf argument for each dataset
    config = []  # each element is a dict, keys = {'in1: [], in2:[], 'gtf': str, 'lb': str}
    for isize, fn in files_for_size:
        config.append({'in1': [fn], 'in2': None, 'gtf': None, 'lb': isize})

    # for different gtf files
    config_case_ctrl_base = {'in1': files_case_ctrl['ctrl'], 'in2': files_case_ctrl['case'], 'gtf': None, 'lb': None}
    gtf_list = ['hg19.ensGene.gtf.gz', 'hg19.knownGene.gtf.gz', 'hg19.ncbiRefSeq.gtf.gz', 'hg19.refGene.gtf.gz']
    gtf_list = [f'{pw_gtf}/{fn}' for fn in gtf_list]
    for gtf in gtf_list:
        config_case_ctrl = config_case_ctrl_base.copy()
        config_case_ctrl['gtf'] = gtf
        gtf_lb = os.path.basename(gtf).split('.')[1]
        config_case_ctrl['lb'] = gtf_lb
        config.append(config_case_ctrl)
    
    # for different sample count by simply repeat the case and ctrl files (create symlink by adding _n suffix), 
    # gtf, use the default one (gtf=None)
    # first , create the symlink to sam_dup folder
    os.makedirs('sam_dup', exist_ok=True)
    for n_sam in [2, 4, 8, 16]:
        n_rep = n_sam // 2 
        config_case_ctrl = config_case_ctrl_base.copy()
        config_case_ctrl['lb'] = f'{n_sam}_sam'
        for case_ctrl, v in files_case_ctrl.items():
            in1in2 = {'ctrl': 'in1', 'case': 'in2'}[case_ctrl]
            fls = []
            for fn in v:
                for i in range(n_rep):
                    i += 1
                    suffix = f'_{i}' if n_rep > 1 else ''
                    fn_new = f'sam_dup/{os.path.basename(fn).replace(".bed", f"{suffix}.bed")}'
                    # logger.info(fn_new)
                    if not os.path.exists(fn_new):
                        os.symlink(fn, fn_new)
                    fls.append(fn_new)
            config_case_ctrl[in1in2] = fls
            
        config.append(config_case_ctrl)
    
    # tmp = json.dumps(config, indent=4)
    # logger.info(f'config = {tmp}')
    return config

def create_sh(pwd, pw_python, pw_nrsa_code, config, n_rep):
    """
    create the script for run both step1 and step2, each one repeat n_rep times
    """
    record = {'step1': {}, 'step2': {}}  # k1 = step1, step2, k2 = lb, v = [], each element is like [fn_sh, fn_log]
    organism = 'hg19'
    for c in config:
        lb = c['lb']
        lb_in_fn = re.sub(r'\s+', '_', lb)
        for step in ['step1', 'step2']:
            tmp = {'step1': 'pause_PROseq.py', 'step2': 'eRNA.py'}
            fn_nrsa_script = f'{pw_nrsa_code}/{tmp[step]}'
            for i_rep in range(n_rep):
                i_rep += 1
                sh = f'{pwd}/sh/{lb_in_fn}_{step}_{i_rep}.sh'
                log = f'{pwd}/sh/{lb_in_fn}_{step}_{i_rep}.log'
                record[step].setdefault(lb, []).append([sh, log])

                with open(sh, 'w') as f:
                    f.write(f'#!/bin/bash\n')
                    f.write(f'cd {pwd}\n')
                    f.write(f'echo "start {step} {lb} {i_rep}"\n')
                    f.write(f'/usr/bin/time -v {pw_python}  {fn_nrsa_script} -org {organism} ')
                    if c['in1']:
                        f.write(f' -in1 {" ".join(c["in1"])} ')
                    if c['in2']:
                        f.write(f' -in2 {" ".join(c["in2"])} ')
                    if c['gtf']:
                        f.write(f' -gtf {c["gtf"]} ')
                    f.write(f'-pwout {pwd}/out/{lb_in_fn}_{step}_{i_rep} ')
                    f.write(f'> {log} 2>&1\n')
                    f.write(f'echo "end {step} {lb} {i_rep}"\n')
    # dump the record to json
    with open(f'{pwd}/bench_settings.json', 'w') as f:
        json.dump(record, f, indent=4)
    
    # build the sh list for run them sequentially, step1 scripts first, then step2, after each script, wait 30s
    sh_list = []
    for step in ['step1', 'step2']:
        for lb, v in record[step].items():
            for fn_sh, fn_log in v:
                sh_list.append(fn_sh)
                os.chmod(fn_sh, 0o755)
                sh_list.append(f'echo "waiting 30s"\nsleep 30\n')
    sh_list.append('echo "all done"')
    with open(f'{pwd}/run_all.sh', 'w') as f:
        f.write('\n'.join(sh_list))
    
    
if __name__ == "__main__":
    if hostname not in ['cqs3', 'bioinfo2']:
        print('This script should be run in cqs3 or bioinfo2')
        sys.exit(1)
    n_rep = 3 # number of repeats for each dataset, to get the average time and memory usage
    if hostname == 'cqs3':
        pwd = '/nobackup/h_vangard_1/chenh19/nrsa/bench'
        pw_python = '/data/cqs/chenh19/project/nrsa_v2/miniconda3/bin/python3.12'
        pw_nrsa_code = '/data/cqs/chenh19/project/nrsa_v2/code'
    else:
        pwd = '/data/nrsa/bench'
        pw_python = ''
    os.chdir(pwd)
    logger.info('prepare config')
    config = prepare_config(pwd)
    logger.info('create sh')
    create_sh(pwd, pw_python, pw_nrsa_code, config, n_rep)
    