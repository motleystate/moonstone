"""Build, check and write the configuration of a reads run, and read the site file."""
import glob
import os
import re
import shutil
import sys

import yaml

from moonstone.validation import Finding

# The menu of hosts whose reads can be removed: name -> NCBI taxonomy ID
HOSTS = {'human': 9606, 'mouse': 10090}
EXTRACT_READS = 'extract_kraken_reads.py'     # from KrakenTools
SITE_FILE_VARIABLE = 'MOONSTONE_SITE'


def site() -> dict:
    """
    Site settings, from the file named by $MOONSTONE_SITE, else ~/.config/moonstone/site.yaml; {} without a file.

    Keys, all optional: ``db`` and ``host_db`` (default Kraken2 databases), ``adapter_dir`` (Trimmomatic
    adapters), ``submit`` (command that submits a Slurm job, default ``sbatch``), ``head_options`` (its options
    for the job that runs Snakemake) and ``profile`` (a Snakemake Slurm profile folder).
    """
    path = os.environ.get(SITE_FILE_VARIABLE) or os.path.expanduser('~/.config/moonstone/site.yaml')
    try:
        with open(path) as handle:
            settings = yaml.safe_load(handle)
    except OSError:
        return {}
    return settings if isinstance(settings, dict) else {}


def build(reads_link_dir: str, output_dir: str, paired: bool, ext: str, options: dict) -> dict:
    """
    The workflow configuration (the keys of workflow/config_all_options.yaml) for one run.

    :param options: ``db``, and optionally ``host_db``, ``remove_host`` (names from :data:`HOSTS`; default all,
        empty: no host removal), ``no_trimming``, ``no_bracken``, ``read_len``, ``bracken_threshold``,
        ``keep_trimmed_reads``, ``keep_kraken_output``, ``min_reads``, ``min_classified_pct``, ``adapter_dir``.
    """
    hosts = options.get('remove_host')
    hosts = list(HOSTS) if hosts is None else list(hosts)
    unknown = [host for host in hosts if host not in HOSTS]
    if unknown:
        raise ValueError(f"Unknown host(s): {', '.join(unknown)}. Available: {', '.join(HOSTS)}.")
    settings = {
        'input_data': {'raw_read_dir': reads_link_dir, 'paired': bool(paired), 'ext': ext},
        'to_perform': {'trimmomatic': not options.get('no_trimming', False), 'host_mapping': bool(hosts),
                       'bracken': not options.get('no_bracken', False)},
        'krakenDB_dir': options.get('db'),
    }
    if hosts:
        settings['krakenDB_host'] = options.get('host_db')
        settings['host_taxid'] = [HOSTS[host] for host in hosts]
    settings.update({
        'keep_trimmed_reads': bool(options.get('keep_trimmed_reads', False)),
        'keep_kraken_output': bool(options.get('keep_kraken_output', False)),
        'output_dir': output_dir,
        'bracken': {'read_len': options.get('read_len') or 150, 'threshold': options.get('bracken_threshold') or 5},
        'report': {'min_reads': 100000 if options.get('min_reads') is None else options['min_reads'],
                   'min_classified_pct': 20 if options.get('min_classified_pct') is None
                   else options['min_classified_pct']},
    })
    if options.get('adapter_dir'):
        settings['trimmomatic'] = {'adapter_dir': options['adapter_dir']}
    return settings


def adapter_dir(settings: dict):
    """The Trimmomatic adapter folder: trimmomatic.adapter_dir, else the conda trimmomatic package's, else None."""
    packaged = glob.glob(os.environ.get('CONDA_PREFIX', sys.prefix) + '/share/trimmomatic*/adapters')
    return (settings.get('trimmomatic') or {}).get('adapter_dir') or (packaged[0] if packaged else None)


def validate(settings: dict) -> list:
    """
    Everything that would stop the workflow before its first step: a missing setting, read files not found or
    not named {sample}_1/{sample}_2, a folder that is not a Kraken2 database, a missing Bracken k-mer file,
    missing adapters, or the host-removal program not installed. The workflow itself calls this function.

    :return: list of :class:`moonstone.validation.Finding`, all errors; empty when the run can start.
    """
    def error(code, message):
        return Finding('error', code, message)

    if 'keep_intermediates' in settings:
        return [error('old_setting', 'keep_intermediates is no longer a setting: use keep_trimmed_reads and '
                                     'keep_kraken_output')]
    input_data, perform = settings.get('input_data') or {}, settings.get('to_perform') or {}
    do_trim, do_host = perform.get('trimmomatic', True), perform.get('host_mapping', True)
    reads_dir, database = input_data.get('raw_read_dir'), settings.get('krakenDB_dir')
    host_database = settings.get('krakenDB_host')
    required = [('input_data: raw_read_dir', reads_dir), ('krakenDB_dir', database)]
    if do_host:
        required.append(('krakenDB_host (needed unless to_perform: host_mapping is False)', host_database))
    missing = [error('missing_setting', 'Missing required setting: ' + name) for name, value in required if not value]
    if missing:
        return missing

    findings = []
    paired, ext = input_data.get('paired', True), input_data.get('ext', 'fastq.gz')
    if not os.path.isdir(reads_dir):
        findings.append(error('not_a_directory', 'input_data.raw_read_dir is not a directory: ' + reads_dir))
    else:
        pattern = re.compile((r'(.+)_(.+)\.' if paired else r'(.+)\.') + re.escape(ext))
        matches = [match for match in map(pattern.fullmatch, os.listdir(reads_dir)) if match]
        mates = sorted({'_' + match.group(2) for match in matches}) if paired else ['']
        if not matches:
            findings.append(error('no_reads', 'No input files found in ' + reads_dir
                                  + " matching the 'paired'/'ext' settings"))
        elif paired and set(mates) != {'_1', '_2'}:
            findings.append(error('mate_names', 'paired: True expects {sample}_1/{sample}_2 files, found read suffixes: '
                                  + str(mates)))
    for name, path in [('krakenDB_dir', database)] + ([('krakenDB_host', host_database)] if do_host else []):
        if not os.path.isdir(path):
            findings.append(error('not_a_directory', name + ' is not a directory: ' + path))
        elif not os.path.exists(os.path.join(path, 'hash.k2d')):
            findings.append(error('not_a_database', name + ' is not a Kraken2 database (no hash.k2d): ' + path))
    read_len = (settings.get('bracken') or {}).get('read_len', 150)
    kmers = os.path.join(database, 'database%smers.kmer_distrib' % read_len)
    if perform.get('bracken', True) and os.path.isdir(database) and not os.path.exists(kmers):
        findings.append(error('no_kmer_file', 'Bracken needs %s; set bracken.read_len to an existing length or generate '
                                              'it with bracken-build -d %s -l %s' % (kmers, database, read_len)))
    if do_trim and not adapter_dir(settings):
        findings.append(error('no_adapters', 'Trimmomatic adapters not found: set trimmomatic: adapter_dir, or install '
                                             'the conda trimmomatic package'))
    if do_host and not shutil.which(EXTRACT_READS):
        findings.append(error('no_krakentools', 'Host-read removal needs %s from KrakenTools on PATH (conda package: '
                                                'krakentools)' % EXTRACT_READS))
    return findings


def write(settings: dict, path: str) -> None:
    """Write the configuration as YAML: the record of a run's settings. An existing file is never overwritten."""
    with open(path, 'x') as handle:
        yaml.safe_dump(settings, handle, sort_keys=False)
