#! /usr/bin/env python3
"""
benchmark for different datasets, run in cqs3
"""
import sys
import os, re
import json
import traceback # use trackback.format_exc() to capture info
home = os.path.expanduser("~")
import numpy as np

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
        # '/workspace/wangj52/PROseq-data/Bcell-GM12004-GSM980646/bams/GSM980646-GM12004_PROseq_sorted_rRNArm-F4q10.bam',
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
        fn_bed = re.sub(r'\.bed$', '.sorted.bed', fn_bed)

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

def create_sh(pwd, config, n_rep, version=None, force=False):
    """
    create the script for run both step1 and step2, each one repeat n_rep times
    """
    if version not in {'v1', 'v2'}:
        logger.error(f'invlid version: {version}')
        sys.exit(1)

    if version == 'v1':
        pw_suffix = '_v1'
        organism = '-m hg19'
        arg_extra = ''
        pw_program = 'perl'
        executable_map = {'step1': 'pause_PROseq.pl', 'step2': 'eRNA.pl'}
        pw_nrsa_code  = '/data/nrsa/v1/bin'
        pwout_arg_k = {'step1': '-o', 'step2': '-w'}
    else:
        pw_suffix = ''
        organism = '-org hg19'
        arg_extra = '-v -sorted '
        pw_program = '/data/cqs/chenh19/project/nrsa_v2/miniconda3/bin/python3.12' if hostname == 'cqs3' else '/data/miniconda/envs/nrsa/bin/python'
        executable_map = {'step1': 'pause_PROseq.py', 'step2': 'eRNA.py'}
        pw_nrsa_code = '/data/cqs/chenh19/project/nrsa_v2/code' if hostname == 'cqs3' else '/data/nrsa/code'
        pwout_arg_k = {'step1': '-pwout', 'step2': '-pwout'}

    # perl /data/nrsa/v1/bin/pause_PROseq.pl -m mm10 -o /data/nrsa/testv1/mm10 -in1 /home/chenh19/data/nrsa/testdata/INF_EV_0hr_1-R1.sorted.bed /home/chenh19/data/nrsa/testdata/INF_EV_0hr_2-R1.sorted.bed -in2 /home/chenh19/data/nrsa/testdata/INF_Cre_0hr_1-R1.sorted.bed /home/chenh19/data/nrsa/testdata/INF_Cre_0hr_2-R1.sorted.bed  &> /data/nrsa/bench/v1.step1.mm10.log;
    # sleep 60;
    # /usr/bin/time -v perl /data/nrsa/v1/bin/eRNA.pl -m mm10 -w  /data/nrsa/testv1/mm10 -in1 /home/chenh19/data/nrsa/testdata/INF_EV_0hr_1-R1.sorted.bed /home/chenh19/data/nrsa/testdata/INF_EV_0hr_2-R1.sorted.bed -in2 /home/chenh19/data/nrsa/testdata/INF_Cre_0hr_1-R1.sorted.bed /home/chenh19/data/nrsa/testdata/INF_Cre_0hr_2-R1.sorted.bed &> /data/nrsa/bench/v1.step2.mm10.log;

    # /usr/bin/time -v perl /data/nrsa/v1/bin/pause_PROseq.pl -m hg19  -in1 sam_dup/4104-AW-1.sorted_1.bed sam_dup/4104-AW-1.sorted_2.bed sam_dup/4104-AW-1.sorted_3.bed sam_dup/4104-AW-1.sorted_4.bed sam_dup/4104-AW-3.sorted_1.bed sam_dup/4104-AW-3.sorted_2.bed sam_dup/4104-AW-3.sorted_3.bed sam_dup/4104-AW-3.sorted_4.bed  -in2 sam_dup/4104-AW-2.sorted_1.bed sam_dup/4104-AW-2.sorted_2.bed sam_dup/4104-AW-2.sorted_3.bed sam_dup/4104-AW-2.sorted_4.bed sam_dup/4104-AW-4.sorted_1.bed sam_dup/4104-AW-4.sorted_2.bed sam_dup/4104-AW-4.sorted_3.bed sam_dup/4104-AW-4.sorted_4.bed -o /data/nrsa/bench/out_v1/8_sam_step1_2 > /data/nrsa/bench/sh_v1/8_sam_step1_2.log 2>&1
    record = {'step1': {}, 'step2': {}}  # k1 = step1, step2, k2 = lb, v = [], each element is like [fn_sh, fn_log]
    
    for c in config:
        lb = c['lb']
        lb_in_fn = re.sub(r'\s+', '_', lb)
        
        if version == 'v1':
            if lb == '16_sam' or c['gtf']:
                continue
        
        
        for step in ['step1', 'step2']:
            fn_nrsa_script = f'{pw_nrsa_code}/{executable_map[step]}'
            for i_rep in range(n_rep):
                i_rep += 1
                sh = f'{pwd}/sh{pw_suffix}/{lb_in_fn}_{step}_{i_rep}.sh'
                log = f'{pwd}/sh{pw_suffix}/{lb_in_fn}_{step}_{i_rep}.log'

                record[step].setdefault(lb, []).append([sh, log])

                with open(sh, 'w') as f:
                    f.write(f'#!/bin/bash\n')
                    f.write(f'cd {pwd}\n')
                    pwout_arg = pwout_arg_k[step]
                    pwout = f'{pwd}/out{pw_suffix}/{lb_in_fn}_step1_{i_rep}'
                    rm_command = f'\nrm -rf {pwout}' if step == 'step1' else ''
                    f.write(f'echo "start {step} {lb} {i_rep}"{rm_command}\ndate\n')
                    f.write(f'/usr/bin/time -v {pw_program}  {fn_nrsa_script} {organism} {arg_extra}')
                    
                    
                    if c['in1']:
                        f.write(f' -in1 {" ".join(c["in1"])} ')
                    if c['in2']:
                        f.write(f' -in2 {" ".join(c["in2"])} ')
                    if c['gtf'] and version == 'v2':
                        f.write(f' -gtf {c["gtf"]} ')
                    f.write(f'{pwout_arg} {pwout} ')
                    f.write(f'> {log} 2>&1\nif [[ $? -ne 0 ]]; then\n\techo "error in {step} {lb} {i_rep}" >>{log}\nelse\n\techo "script finished without error" >>{log}\nfi\n')
                    f.write(f'echo "end {step} {lb} {i_rep}"\ndate\n')
    # dump the record to json

    
    # build the sh list for run them sequentially, step1 scripts first, then step2, after each script, wait 30s
    sh_list = []
    sleep_time = 10
    ct = {'already_done': 0, 'failed': 0, 'new': 0}
    for step in ['step1', 'step2']:
        for lb, v in record[step].items():
            for fn_sh, fn_log in v:
                if not force and os.path.exists(fn_log):
                    # check if there are "script finished without error" in the log file, if yes, skip
                    with open(fn_log) as f:
                        if 'script finished without error' in f.read():
                            sh_list.append(f'echo "already done - skip {fn_sh}"')
                            ct['already_done'] += 1
                            continue
                        else:
                            logger.warning(f'log file exist but not finished: {fn_log}')
                            ct['failed'] += 1
                else:
                    ct['new'] += 1
                
                
                sh_list.append(fn_sh)
                os.chmod(fn_sh, 0o755)
                sh_list.append(f'echo "waiting {sleep_time}s"\nsleep {sleep_time}\n')
    sh_list.append('echo "all done"')
    
    logger.info(f'version = {version} - {ct}')
    
    with open(f'{pwd}/run_all_{version}.sh', 'w') as f:
        f.write('\n'.join(sh_list))
    return record

def parse_nrsa_log(fn_log):
    ts_count = 'NA'
    mem = None
    time = None
    run_success = False
    with open(fn_log) as f:
        for i in f:
            if 'initial gtf dict =' in i:
                m = re.search(r'current ts count = (\d+)', i)
                ts_count = int(m.group(1))
            if 'Elapsed (wall clock) time' in i:
                # line is like  Elapsed (wall clock) time (h:mm:ss or m:ss): 0:30.28 , convert the time to seconds, can also include the hour 
                time = i.strip().rsplit(' ', 1)[-1]
                time_split = [float(_) for _ in time.split(':')]
                # fill the time to 3 elements, if only 1 element, fill the first 2 with 0
                time_split = [0] * (3 - len(time_split)) + time_split
                time = time_split[-1] + time_split[-2] * 60 + time_split[-3] * 3600
            elif 'Maximum resident set size (kbytes):' in i:
                mem = int(i.strip().rsplit(' ', 1)[-1]) / 1024 / 1024  # convert to GB
            elif 'script finished without error' in i:
                run_success = True
    # if ts_count == 'NA':
    #     logger.error(f'ts_count not found: {fn_log}')
    if mem is None or time is None:
        logger.error(f'error in parsing log file: {fn_log}')
        return ts_count, None, None
    
    if not run_success:
        logger.error(f'error in running script: {fn_log}')
        return ts_count, None, None
    
    
    return ts_count, time, mem

def extract_number(s):
    return float(re.search(r'[\d\.]+', s).group())

def plot_benchmark_results(df, fn_figure):
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    import seaborn as sns
    
    category_map = {'size': 'Reads count', 'sam_number': 'Sample Number', 'gtf': 'GTF'}
    
    fig = plt.figure(figsize=(15, 20))
    text_row_height = 0.05
    plot_row_height = 1.1
    gs = GridSpec(6, 2, height_ratios=[text_row_height, plot_row_height, text_row_height, plot_row_height, text_row_height, plot_row_height], figure=fig)
    # fig, axes = plt.subplots(3, 2, figsize=(15, 20))
    # plt.subplots_adjust(hspace=0.5, top=0.95)
    tick_label_size=18
    tick_label_size_gtf_row = 15
    axis_title_size=21
    title_size=25
    bar_label_padding = 8
    
    from matplotlib.container import ErrorbarContainer
    
    for i, category in enumerate(df['category'].unique()):
        data = df[df['category'] == category]
        # Calculate mean and std for time and memory
        grouped = data.groupby('lb').agg({'time': ['mean', 'std'], 'mem': ['mean', 'std']})
        
        grouped = grouped.sort_index(key=lambda x: x.map(extract_number))

        if category == 'gtf':
            x_labels = [_.replace(' ', '\n') for _ in grouped.index]
        else:
            x_labels = grouped.index

        ax_text = fig.add_subplot(gs[i * 2, :])
        ax_text.text(0.5, -0.4, category_map[category], 
                    ha='center', va='center', fontsize=title_size, fontweight='bold',
                    bbox=dict(facecolor='none', pad=0, edgecolor='none'))
        ax_text.axis('off')  # Hide axes for text box
        # for spine in ax_text.spines.values():
        #     spine.set_visible(True)
        #     spine.set_color('red')

        # Time plot
        # ax = axes[i, 0]
        ax = fig.add_subplot(gs[i * 2 + 1, 0])
        # x = range(len(grouped.index))
        # ax.bar(x, grouped[('time', 'mean')], yerr=grouped[('time', 'std')], capsize=5)
        # use seaborn to plot barplot with hue
        # ax.set_xticks(x)
        # ax.set_xticklabels(grouped.index)
        
        sns.barplot(x='lb', y='time', hue='version', data=data, 
                capsize=0.1, errwidth=1, errorbar='sd', ax=ax, order=grouped.index)
        ax.set_ylabel('Time (s)', fontsize=axis_title_size)
        ax.tick_params(axis='y', labelsize=tick_label_size)
        ax.tick_params(axis='x', labelsize=tick_label_size_gtf_row if category == 'gtf' else tick_label_size)
        ax.set_xticklabels(x_labels)
        ax.set_xlabel('')
        
        error_bars = ax.lines
        
        idx_bar = 0
        for container in ax.containers:
            labels = ax.bar_label(container, padding=7, fmt='%.0f')
            plt.setp(labels, fontsize=tick_label_size - 1)

            # n_labels = len(labels)
            # for label in labels:
            #     error_bar_top = error_bars[idx_bar].get_ydata()
            #     idx_bar += 1
            #     if not np.isnan(error_bar_top[1]):
            #         label.set_y(error_bar_top[1])
            #         logger.info(error_bar_top[1])

        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax * 1.05)
        
        
        # Memory plot
        # ax = axes[i, 1]
        ax = fig.add_subplot(gs[i * 2 + 1, 1])
        # ax.bar(x, grouped[('mem', 'mean')], yerr=grouped[('mem', 'std')], capsize=5)
        sns.barplot(x='lb', y='mem', hue='version', data=data, 
                capsize=0.1, errwidth=1, errorbar='sd', ax=ax, order=grouped.index)
        ax.set_ylabel('Memory (GB)', fontsize=axis_title_size)
        ax.tick_params(axis='y', labelsize=tick_label_size)
        ax.tick_params(axis='x', labelsize=tick_label_size_gtf_row if category == 'gtf' else tick_label_size)
        ax.set_xticklabels(x_labels)
        ax.set_xlabel('')
        for container in ax.containers:
            labels = ax.bar_label(container, padding=bar_label_padding, fmt='%.2f')
            plt.setp(labels, fontsize=tick_label_size - 1)
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax * 1.05)

        # ax.set_xticks(x)
        # ax.set_xticklabels(grouped.index)
        
        # for spine in ax.spines.values():
        #     spine.set_visible(True)
        #     spine.set_color('red')

        # # Add category label as row title
        # fig.text(0.01, 0.77 - i*0.31, category_map[category], rotation=90, 
        #          ha='left', va='center', fontsize=title_size, fontweight='bold')
    
    # # Set column titles
    # fig.text(0.3, 0.95, 'Time (s)', ha='center', va='bottom', fontsize=title_size, fontweight='bold')
    # fig.text(0.75, 0.95, 'Memory (GB)', ha='center', va='bottom', fontsize=title_size, fontweight='bold')
    
    # plt.subplots_adjust(left=0, bottom=0.08, right=1, top=0.92, wspace=0.01, hspace=0.01)
    plt.tight_layout(rect=[0.03, 0, 1, 1], h_pad=0)  # Adjust layout to accommodate row titles
    # plt.tight_layout()  # Adjust layout to accommodate row titles
    plt.savefig(fn_figure, dpi=300, bbox_inches='tight')

def collect_bench_res(script_settings):
    import pandas as pd
    res_l = []
    n_rep = 3
    
    lb_convert = {}
    
    # fn_bam_size = f'{pwd}/bam_file_size.json'
    # with open(fn_bam_size, 'r') as f:
    #     bam_size = json.load(f)
    #     for k, v in bam_size.items():
    #         if not k.endswith('.sorted.bed'):
    #             continue
    #         size, reads = v['size'], v['reads']
    #         if size not in lb_convert:
    #             lb_convert[size] = reads
    #         else:
    #             logger.error(f'duplicate size: {size}')

    not_ready = []
    ts_count_d = {}
    res_step1_step2 = {}
    for version in ['v1', 'v2']:
        config = script_settings[version]
        
        # logger.warning(f'modify here, use version2 logs to mimic version1')
        # config = script_settings['v2']
        
        for step in ['step1', 'step2']:
            v1 = config[step]
            for lb, v2 in v1.items():
                # if lb ends with GB, category = size, if ends with _sam, category = sam_number, else category = gtf
                category = 'size' if lb.endswith('GB') else 'sam_number' if lb.endswith('_sam') else 'gtf'
                for i_rep, (fn_sh, fn_log) in enumerate(v2):
                    # if fn_log.endswith(f'16_sam_step1_3.log'):
                    #     fn_log = fn_log.replace('step1_3.log', 'step1_3.new.log')
                    
                    if not os.path.exists(fn_log):
                        not_ready.append(fn_log)
                        continue
                    
                    ts_count, time, mem = parse_nrsa_log(fn_log)
                    if ts_count != 'NA':
                        ts_count_d[lb] = ts_count
                    else:
                        ts_count = ts_count_d.get(lb, 'NA')
                    
                    if time is None or mem is None:
                        not_ready.append(fn_log)
                        continue
                    
                    if version == 'v1' and category == 'gtf':
                        continue
                    if version == 'v1' and lb == '16_sam':
                        continue
                    
                    # create the new label
                    lb_new = lb_convert.get(lb, lb)
                    if category == 'gtf':
                        ts_count_new = f'{ts_count/1000:.1f}K'
                        lb_new = f'{lb} {ts_count_new}'
                    if category == 'sam_number':
                        lb_new = f'n={lb.split("_")[0]}'
                    res_l.append([version, step, category, lb_new, i_rep, time, mem])
                    
                    k_step1_step2 = (version, category, lb_new, i_rep)
                    if k_step1_step2 not in res_step1_step2:
                        res_step1_step2[k_step1_step2] = [version, 'step1_step2', category, lb_new, i_rep, 0, 0]
                    res_step1_step2[k_step1_step2][-2] += time
                    res_step1_step2[k_step1_step2][-1] = max(res_step1_step2[k_step1_step2][-1], mem)
        
    for v in res_step1_step2.values():
        res_l.append(v)
    
    if not_ready:
        tmp = '\n\t'.join(not_ready)
        logger.error(f'files not ready n = {len(not_ready)}:\n\t{tmp}')

    # write to csv
    with open('bench_res.csv', 'w') as f:
        f.write('version,step,category,lb,i_rep,time,mem\n')
        for r in res_l:
            f.write(','.join(map(str, r)) + '\n')

    # convert to pandas and plot
    df = pd.DataFrame(res_l, columns=['version', 'step', 'category', 'lb', 'i_rep', 'time', 'mem'])
    for k in ['step1', 'step2', 'step1_step2']:
        df_k = df[df['step'] == k]
        
        fn_figure = f'bench_res_{k}.png'
        plot_benchmark_results(df_k, fn_figure)
        logger.info(f'figure saved to {fn_figure}')

    
if __name__ == "__main__":
    import argparse as arg
    from argparse import RawTextHelpFormatter
    ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    ps.add_argument('action', help="""action to run, can be script/build, collect/plot""", choices=['script', 'build', 'collect', 'plot'])
    ps.add_argument('-force', '-f', help="""force build the script, ignore the existing log files""", action='store_true')
    args = ps.parse_args()

    force = args.force
    if hostname not in ['cqs3', 'bioinfo2']:
        print('This script should be run in cqs3 or bioinfo2')
        sys.exit(1)
    n_rep = 3 # number of repeats for each dataset, to get the average time and memory usage
    if hostname == 'cqs3':
        pwd = '/nobackup/h_vangard_1/chenh19/nrsa/bench'
    else:
        pwd = '/data/nrsa/bench'
    os.chdir(pwd)
    
    
    # create folder
    for ipw in ['sh', 'sh_v1', 'out', 'out_v1']:
        os.makedirs(f'{pwd}/{ipw}', exist_ok=True)
    
    action = args.action
    # collapse actions, script = build, collect = plot
    if action == 'script':
        action = 'build'
    elif action == 'collect':
        action = 'plot'

    logger.info('prepare config')
    config = prepare_config(pwd)
    
    
    if action == 'build':
        logger.info('create sh')
        script_settings = {}
        for version in ['v1', 'v2']:
            script_settings[version] = create_sh(pwd, config, n_rep, version=version, force=force)
        with open(f'{pwd}/bench_settings.json', 'w') as f:
            json.dump(script_settings, f, indent=4)

    elif action == 'plot':
        fn_json = f'{pwd}/bench_settings.json'
        with open(fn_json, 'r') as f:
            script_settings = json.load(f)
        collect_bench_res(script_settings)
