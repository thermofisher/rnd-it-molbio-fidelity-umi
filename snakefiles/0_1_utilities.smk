def get_sample_files_script(sample, debug=True):

    input_dir = config['MAIN__INPUT_DIR']

    r1_files = []
    r1_stems = []
    r2_files = []
    r2_stems = []
    r1_testlist = []
    r2_testlist = []
    sfq = False

    wlogger.info(f'SEARCHING FOR FASTQ:\tDIR\t{config["MAIN__INPUT_DIR"]}\tSAMPLE\t{sample}')

    try:
        wlogger.info(f'SEARCHING FOR FASTQ:\tChecking if `slimfastq` is configured.')
        subprocess.call([f'{CURRENT_PATH}/dependencies/slimfastq/slimfastq', '-h'], stdout = subprocess.DEVNULL)
    except:
        wlogger.error(f'SEARCHING FOR FASTQ:\t`slimfastq` was not found, Run this cmd to install.')
        wlogger.error(f'git submodule update --init --recursive; cd {CURRENT_PATH}/dependencies/slimfastq; git checkout v2.01 && make; cd ../../;')
        sys.exit(1)

    if debug:
        wlogger.info(f'SEARCHING FOR FASTQ:\tLANE REGEX:\t{config["MAIN__LANE_RE"]}')
    for root, dirnames, filenames in os.walk(input_dir):
        for f in filenames:
            if ( f.endswith('.fastq.gz') or f.endswith('.fq.gz') or f.endswith('.sfq') ) and f.find(sample) > -1:

                if f.find('.sfq') > -1:
                    sfq = True
                if debug:
                    wlogger.info(f'File found for {sample}: {f} ')

                if f.find('_R1_') > -1 or f.endswith('_1.fq.gz') or f.endswith('_1.sfq') or f.endswith('_1.fastq.gz'):
                    if debug:
                        wlogger.info(f'SEARCHING FOR FASTQ:\tDetected R1 file:\t{f}\tin\t{root}')
                    r1_files.append(root+os.sep+f)
                    if config['MAIN__LANE_RE'] != "":
                        lane_splits = re.split(config['MAIN__LANE_RE'], f)
                        if lane_splits[0].find(sample) > -1:
                            stem = lane_splits[0]
                        else:
                            stem = lane_splits[1]
                    else:
                        stem = sample
                    r1_testlist.append(root + os.sep + stem)
                    r1_stems.append(stem)

                elif f.find('_R2_') > -1 or f.endswith('_2.fq.gz') or f.endswith('_2.sfq') or f.endswith('_2.fastq.gz'):
                    if debug:
                        wlogger.info(f'SEARCHING FOR FASTQ:\tDetected R2 file:\t{f}\tin\t{root}')
                    r2_files.append(root + os.sep + f)
                    if config['MAIN__LANE_RE'] != "":
                        lane_splits = re.split(config['MAIN__LANE_RE'], f)
                        if lane_splits[0].find(sample) > -1:
                            stem = lane_splits[0]
                        else:
                            stem = lane_splits[1]
                    else:
                        stem = sample
                    r2_testlist.append(root + os.sep + stem)
                    r2_stems.append(stem)

    all_stems = []
    all_stems.extend(r1_stems)
    all_stems.extend(r2_stems)

    if len(list(set(all_stems))) > 1:
        wlogger.warning(f'SEARCHING FOR FASTQ:\tSAMPLE\t{sample}\tis matching to more than one sample: {all_stems}')

    r1_testlist.sort()
    r2_testlist.sort()
    if str(r1_testlist) != str(r2_testlist):
        wlogger.error(f'SEARCHING FOR FASTQ:\tPaired-end mode, yet found R1 and R2 lists does not match!')
        wlogger.error(f'SEARCHING FOR FASTQ:\tR1 starts\t{r1_testlist}')
        wlogger.error(f'SEARCHING FOR FASTQ:\tR2 starts\t{r2_testlist}')
        raise Exception(f'Key {sample} is not properly differentiating, pairs are invalid\nCheck lists of R1: {str(r1_testlist)}\nand R2: {str(r2_testlist)}')

    if len(r1_files) == 0:
        if not os.path.exists(input_dir):
            wlogger.error(f'SEARCHING FOR FASTQ:\tFileExistsError\t{input_dir}')
        else:
            wlogger.error(f'SEARCHING FOR FASTQ:\tFileNotFoundError\t{input_dir}\t{sample}')
            list_dir = os.listdir(input_dir)
            wlogger.error(f'SEARCHING FOR FASTQ:\tListDir\t{list_dir}')

        raise Exception(f'No matching R1 files for {sample} in {input_dir}')

    r1_files.sort()
    r2_files.sort()
    r1ft = os.path.join(tmp+f'/{sample}R1_input_list.sh')
    r2ft = os.path.join(tmp+f'/{sample}R2_input_list.sh')
    f1 = open(r1ft, 'w')

    def shstring(l, sfq):
        slimfq = f'{CURRENT_PATH}/dependencies/slimfastq/slimfastq '
        out = ''
        for i in l:
            i = i.replace(' ', '\\ ')
            if i.endswith('.sfq') and sfq:
                out = f'{out}{slimfq}{i} | gzip -f; \n'
            else:
                out = f'{out}cat {i}; \n'
        return out

    f1.write(shstring(r1_files, sfq))
    f1.close()
    os.system(f'chmod a+x {r1ft}')

    f2 = open(r2ft, 'w')
    out = shstring(r2_files, sfq)
    f2.write(out)
    f2.close()
    os.system(f'chmod a+x {r2ft}')


def get_threads(mem):
    if mem > max_mem:
        mem = max_mem

    threads = int(max_threads / (max_mem / mem))

    return threads


def update_config_values(updatable_dict, provided_dict):
    """To set default parameters if not provided in configuration."""
    def handle_exception():
        wlogger.error(f'It seems that you "{USER}" messed sth up in your configuration file!')
        wlogger.error(f'KeyError: {k} not found in default workflow configuration dictionary.')
        wlogger.error(f'Make sure that you are using a correct yaml template for this. You should probably check a `default.yaml` file for available options.')
        sys.exit()

    for k, v in provided_dict.items():
        if isinstance(v, collections.MutableMapping):
            try:
                update_config_values(updatable_dict[k], provided_dict[k])
            except KeyError as e:
                handle_exception()
        else:
            try:
                if updatable_dict[k] != provided_dict[k]:
                    wlogger.info(f'Updating default configuration parameters: [{k}:{updatable_dict[k]} -> {k}:{provided_dict[k]}]')
                    updatable_dict[k] = provided_dict[k]
            except KeyError as e:
                handle_exception()


def init_log():
    try:
        with open(LOG_CONFIG, 'rt') as f:
            log_config = yaml.safe_load(f.read())
        for handler in log_config['handlers']:
            if log_config['handlers'][handler]['class'] == 'logging.FileHandler':
                fl = log_config['handlers'][handler]['filename']
                log_config['handlers'][handler]['filename'] = LOG_PATH[fl]
                log_dir = LOG_PATH[fl].split('/')[:-1]
                log_dir = '/'.join(log_dir)
                if not os.path.isdir(log_dir):
                    os.makedirs(log_dir)

        logging.config.dictConfig(log_config)

    except Exception as e:
        print(f'{time.asctime()}\tERROR\tLOGS creation failed: {e}')
        sys.exit()


def check_if_reference():

    if config['MAIN__REFERENCE'] == None or config['MAIN__REFERENCE'] == '':
        return False
    else:
        return True


def assert_config(config, must_params):

    wlogger.info(f'ASSERT CONFIG:\tSearching if all mandatory parameters are provided in configuration file.')
    for i in must_params:
        if config[i] is None:
            wlogger.error(f'ASSERT CONFIG:\tAsserting configuration must have parameters failed: [{i}]')
            raise TypeError(f'ASSERT CONFIG:\tFaund NoneType when this paramer is mandatory: {i}')


def min_threads(config_threads, threads):

    try:
        threads = int(threads)
    except ValueError as e :
        wlogger.error(f'CONFIGURATION:\tParsing threads form main command raised an error:\t[{e}]')
        return

    for k, v in config_threads.items():
        if isinstance(v, collections.MutableMapping):
            min_threads(config_threads[k], threads)
        else:
            if config_threads[k] > threads:
                wlogger.info(f'Updating configuration threads: [{k}:{config_threads[k]} -> {k}:{threads}]')
                config_threads[k] = threads


def check_index_file(index_file):
    if not os.path.exists(index_file):
        wlogger.error(f'Index file [{index_file}] does not exits! Please run get_index rule to generate index files!')
        sys.exit(1)
    else:
        f_info = os.stat(index_file)
        wlogger.info(f'Index file found. [{index_file}]')
        wlogger.info(f'Index file (mode, ino, dev, nlink, uid, gid, size, atime, mtime, ctime]:')
        wlogger.info(f'{f_info}')


def flatten(d, parent_key='', sep='__'):
    items = []
    for k, v in d.items():
        valid = re.compile(r"[a-zA-Z_]\w*$")
        if not valid.match(k):
            print(f'Key not mathcing Snakemake regex validation. [{k}]')
            sys.exit(1)
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, collections.MutableMapping):
            items.extend(flatten(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))

    return dict(items)


def ref_size(reference):
    Gb = os.path.getsize(reference)*0.000000001
    wlogger.info(f'REFERENCE SIZE:\t{reference}\t{Gb}Gb')
    if Gb > 0.5:
        size = 'big'
    else:
        size = 'small'
    wlogger.info(f'REFERENCE SIZE:\tset to\t{size}')

    return size


def get_trimm_samples(wildcards):
    if config['UMI_TOOLS__run']:
        return (rules.umi_tools_extract.output[0], rules.umi_tools_extract.output[1])
    else:
        return (rules.get_fastq1.output[0], rules.get_fastq2.output[0])

def get_force_trimm_samples(wildcards):
    if config['TRIMMING__program'] == 'bbduk':
        return (tmp + "/{stem}_R1_001_bbduk.fastq.gz", tmp + "/{stem}_R2_001_bbduk.fastq.gz")
    else:
        return (tmp + "/{stem}_R1_001_adapt.fastq.gz", tmp + "/{stem}_R2_001_adapt.fastq.gz")


def get_processed_samples(wildcards):
    if config['TRIMMING__program'] == 'bbduk' and config["TRIMMING__bbduk__FORCE_TRIM_OPTIONS"]:
        return (tmp + "/{stem}_R1_001_2.fastq.gz", tmp + "/{stem}_R2_001_2.fastq.gz")
    else:
        return get_force_trimm_samples(wildcards)


def get_multiqc_pre(wildcards):
    if config["TRIMMING__bbduk__FORCE_TRIM_OPTIONS"]:
        return expand(multiqc_raw + "/fastqc_report_htmls_zips/{stem}_R{R}_001_2_fastqc.zip", stem=stems, R=Rlist)
    elif config['TRIMMING__program'] == 'bbduk':
        return expand(multiqc_raw + "/fastqc_report_htmls_zips/{stem}_R{R}_001_bbduk_fastqc.zip", stem=stems, R=Rlist)
    else:
        return expand(multiqc_raw + "/fastqc_report_htmls_zips/{stem}_R{R}_001_adapt_fastqc.zip", stem=stems, R=Rlist)