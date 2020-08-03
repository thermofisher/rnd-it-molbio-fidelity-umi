#------------------System configuration file ------------------------#
USER = getpass.getuser()
ARGS = sys.argv
WORK_DIR = os.getcwd()

CURRENT_PATH = "."

if not config:
    print(f'You did not provided mandatory configuration file for this workflow.')
    sys.exit(1)
else:
    config = flatten(config)

# API patameters retrieved from config
import inspect
API = False
for i in inspect.stack():
    if i.function == 'run_snakemake_api':
        API = True
        break
for i in inspect.stack():
    if i.function == 'snakemake':
        fn_args, args_paramname, kwargs_paramname, val_dict = inspect.getargvalues(i.frame)
        if API:
            ARGS.append('-j')
            ARGS.append(val_dict['cores'])
            ARGS.append('--cluster')
            ARGS.append(f'{val_dict["cluster"]}')
        CURRENT_PATH = os.path.dirname(val_dict['snakefile'])
        API == True
        break

for k in config:
    # expecting additional_params sring which if None usualy breaks pipeline
    if 'add' in k:
        if config[k] is None:
            config[k] = ''

# Logger setup
if config['MAIN__OUT_DIR'][-1] == '/':
    out = config['MAIN__OUT_DIR'][0:-1]
else:
    out = config['MAIN__OUT_DIR']

logs = out + '/LOGS'
LOG_PATH = {
    'WORKFLOW_LOG_FILE': logs + '/SNAKEMAKE/workflow.log'
    }
LOG_CONFIG = f'{CURRENT_PATH}/LOG_CONFIG.yaml'
init_log()
wlogger = logging.getLogger('custom_workflow')

# Define temporary and output directories.
if config['MAIN__TMP_DIR'][-1] == '/':
    tmp = config['MAIN__TMP_DIR'][0:-1]
else:
    tmp = config['MAIN__TMP_DIR']

if config['MAIN__BENCHMARKING'] is not None:
    if config['MAIN__BENCHMARKING'][-1] == '/':
        bench = config['MAIN__BENCHMARKING'][0:-1]
    else:
        bench = config['MAIN__BENCHMARKING']
else:
    bench = 'default'

bench = 'benchmarks/' + bench
trim_dir = out + '/TRIMMING'
multiqc_dir = out + '/MULTIQC'
multiqc_raw = out + '/MULTIQC_RAW'
aligned_dir = tmp + '/ALIGNED'

THREADING_CONFIG_PATH = f'{CURRENT_PATH}/snakefiles/threading.yaml'
DEFAULT_CONFIG_PATH = f'{CURRENT_PATH}/default.yaml'

MUST_PARAMETERS = [
    'MAIN__INPUT_DIR',
    'MAIN__SAMPLES_ID',
    'MAIN__TMP_DIR',
    'MAIN__OUT_DIR',
    'MAIN__LANE_RE',
    'MAIN__REFERENCE',
    ]

is_ref = check_if_reference()
wlogger.info(f'User: {USER}')
wlogger.info(f'Wokrking directory: {WORK_DIR}')
wlogger.info(f'CONFIGURATION:\tUpdating default configuration parameters with provided ones.')
assert_config(config, MUST_PARAMETERS)

# Load threading options
wlogger.info(f'CONFIGURATION:\tCONFIGURATION:\tLoading threading confifuration file:\t[{THREADING_CONFIG_PATH}]')
with open(THREADING_CONFIG_PATH) as f:
    try:
        config_threads = yaml.load(f, Loader=yaml.FullLoader)
    except FileNotFoundError as e:
        wlogger.error(f'CONFIGURATION:\tLoading threading confifuration file\t[{THREADING_CONFIG_PATH}] raised an error:\t[{e}] ')
        sys.exit(1)

wlogger.info(f'CONFIGURATION:\tLoading default workflow confifuration file:\t[{DEFAULT_CONFIG_PATH}]')
with open(DEFAULT_CONFIG_PATH) as f:
    try:
        default_config = yaml.load(f, Loader=yaml.FullLoader)
    except FileNotFoundError as e:
        wlogger.error(f'CONFIGURATION:\tLoading default confifuration file\t[{THREADING_CONFIG_PATH}] raised an error:\t[{e}] ')
        sys.exit(1)

provided_config = config.copy()
config = flatten(default_config)
update_config_values(config, provided_config)


if '-j' in ARGS:
    i = ARGS.index('-j')
    min_threads(config_threads, ARGS[i+1])
elif '--jobs' in ARGS:
    i = ARGS.index('--jobs')
    min_threads(config_threads, ARGS[i+1])
elif '--cores' in ARGS:
    i = ARGS.index('--cores')
    min_threads(config_threads, ARGS[i+1])
else:
    min_threads(config_threads, 1)

config['reference_dna_dir'] = os.path.dirname(config['MAIN__REFERENCE'])[0]
config['main_ref_path'] = os.path.splitext(config['MAIN__REFERENCE'])[0]
config['reference_name'] = os.path.basename(config["MAIN__REFERENCE"]).split(".")[0]

max_threads = config['MACHINE__threads']
max_mem = config['MACHINE__memory']

# to avoid two executions on one cluster, since shared memory is used by default.

# Samples
if len(config['MAIN__SAMPLES_ID']) == 0:
    wlogger.info(f'CONFIGURATION:\tSamples were not provided.')
    sys.exit(1)
else:
    stems = config['MAIN__SAMPLES_ID'].split()

wlogger.info(f'CONFIGURATION:\tSamples\t{stems}')
Rlist = ['1', '2']
is_pe = True
se_f = ""
# Trimming
adapters_bbduk = tmp + '/bbduk_adapters.fa'
adapters_adapt = config['TRIMMING__adapters']

# Memory calculations, constants can be changed
if is_ref:
    SIZE = ref_size(config['MAIN__REFERENCE'])
else:
    SIZE = 'big'

if SIZE == 'big':
    mem_java = 30
    # threading for ref loading jl scripts.
    config['MACHINE__threads_julia'] = get_threads(40)
elif SIZE == 'small':
    mem_java = 6
    config['MACHINE__threads_julia'] = get_threads(8)
wlogger.info(f'CONFIGURATION:\tSetting mem_java to {mem_java} Gb.')
wlogger.info(f'CONFIGURATION:\tSetting threads_julia to {config["MACHINE__threads_julia"]}.')

java_rule_threads = get_threads(mem_java)
wlogger.info(f'CONFIGURATION:\tSetting java_rule_threads to {java_rule_threads}.')
umi_fidelity_output = [
    tmp + '/{stem}_fidelity.csv',
    tmp + '/{stem}_fidelity_distributions.csv',
    tmp + '/{stem}_fidelity_size_distribution.csv',
    tmp + '/{stem}_fidelity_error_distributions_per_read.csv',
    tmp + '/{stem}_fidelity_transitversions.csv',
    tmp + '/{stem}_fidelity_transvers_types.csv',
]

if config['FIDELITY__report_mutations']:
    umi_fidelity_output += [
        tmp + '/{stem}_fidelity_passing_mutations.csv',
        tmp + '/{stem}_fidelity_passing_single_mutations.csv',
    ]
if config['FIDELITY__genes_csv']:
    umi_fidelity_output += [
        tmp + '/{stem}_fidelity_passing_aa_mutations.csv',
        tmp + '/{stem}_fidelity_passing_frameshifts.csv',
        tmp + '/{stem}_fidelity_passing_single_aa_mutations.csv',
    ]

umi_fidelity_params = (
    f'-q {config["FIDELITY__alignment_quality"]} '
    f'-Q {config["FIDELITY__base_quality"]} '
    f'-c {config["FIDELITY__cluster_size"]} '
    f'-m {config["FIDELITY__max_mutations_per_read"]} '
    f'-e {config["FIDELITY__referece_mutation"]} '
    f'-s {config["FIDELITY__umi_split_string"]} '
    f'-p {config["FIDELITY__umi_split_position"]} '
    f'-r {config["MAIN__REFERENCE"]} '
    f'--use-grouped-reads '
    f'-L {config["FIDELITY__min_cluster_counts"]} '
    f'-M {config["FIDELITY__min_mutation_coverage"]}'
)

if config['FIDELITY__report_mutations']:
    umi_fidelity_params = f'{umi_fidelity_params} -R '
if config['FIDELITY__genes_csv']:
    umi_fidelity_params = f'{umi_fidelity_params} -G {config["FIDELITY__genes_csv"]} '
    umi_fidelity_params = f'{umi_fidelity_params} -C {config["FIDELITY__genetic_code"]} '
    umi_fidelity_params = f'{umi_fidelity_params} -T'
if config['FIDELITY__additional_params']:
    umi_fidelity_params = f'{umi_fidelity_params} {config["FIDELITY__additional_params"]}'
if config['UMI_TOOLS__run'] == False:
    calculate_fidelity_sort_by_name_input = [
        tmp + '/{stem}_fixed_mates.bam'
    ]
else:
    calculate_fidelity_sort_by_name_input = [
        tmp + '/{stem}_grouped.bam'
    ]

index_file = ""
index_file = config["MAIN__REFERENCE"] + ".bwt"
check_index_file(index_file)
