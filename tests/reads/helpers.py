"""Helpers for the reads tests: dummy inputs, and fake versions of the pipeline's tools.

The fake tools only create the files the real ones would, so the whole workflow can run in seconds without
Kraken2, Bracken, Trimmomatic or a database.
"""
import importlib.util
import os
import stat
import subprocess
import sys

import yaml

HAS_SNAKEMAKE = importlib.util.find_spec('snakemake') is not None
LINEAGE = ('k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; '
           'g__Escherichia; s__coli')

FAKE_TOOL = r'''#!{python}
import os, shutil, sys
name, a = os.path.basename(sys.argv[0]), sys.argv[1:]
def value(flag): return a[a.index(flag) + 1]
def touch(path): open(path, 'w').close()
if name == 'trimmomatic':
    if a[:1] == ['-version']:
        print('0.fake'); sys.exit()
    files = a[3:-1]                              # mode -threads N <files> ILLUMINACLIP:...
    for path in files[2 if a[0] == 'PE' else 1:]:
        touch(path)
elif name == 'kraken2':
    if '--version' in a:
        print('Kraken version 2.fake'); sys.exit()
    with open(value('--report'), 'w') as f:      # percent, clade reads, direct reads, rank, taxid, name
        f.write('10.00\t100\t100\tU\t0\tunclassified\n90.00\t900\t0\tR\t1\troot\n'
                '30.00\t300\t300\tS\t9606\tHomo sapiens\n60.00\t600\t600\tS\t562\tEscherichia coli\n')
    touch(value('--output'))
elif name == 'extract_kraken_reads.py':
    touch(value('-o'))
    if '-o2' in a:
        touch(value('-o2'))
elif name == 'bracken':
    if a[:1] == ['-v']:
        print('Bracken v3.fake'); sys.exit()
    shutil.copy(value('-i'), value('-w')); touch(value('-o'))
elif name == 'kraken-biom':                      # columns are named after the report files, minus the extension
    reports = a[:a.index('--max')]
    open(value('-o'), 'w').write('\n'.join(os.path.basename(r).rsplit('.', 1)[0] for r in reports))
elif name == 'biom':
    samples = open(value('-i')).read().split()
    with open(value('-o'), 'w') as f:
        f.write('# Constructed from biom file\n#OTU ID\t' + '\t'.join(samples) + '\ttaxonomy\n')
        f.write('562\t' + '\t'.join('600.0' for _ in samples) + '\t{lineage}\n')
'''
TOOLS = ('trimmomatic', 'kraken2', 'extract_kraken_reads.py', 'bracken', 'kraken-biom', 'biom')


def fake_bin(folder, without=()):
    """Write the fake tools in ``folder`` and return a PATH that finds them first."""
    os.makedirs(folder, exist_ok=True)
    for tool in TOOLS:
        if tool not in without:
            path = os.path.join(folder, tool)
            with open(path, 'w') as handle:
                handle.write(FAKE_TOOL.replace('{python}', sys.executable).replace('{lineage}', LINEAGE))
            os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR)
    return folder + os.pathsep + '/usr/bin' + os.pathsep + '/bin'


def sparse_file(path, megabytes):
    with open(path, 'wb') as handle:
        handle.truncate(megabytes * 1024 * 1024)


def make_inputs(tmp, paired=True, samples=('S1', 'S2')):
    """Dummy reads and fake databases (classification 3000 MiB, host 1000 MiB); returns a workflow config."""
    reads = os.path.join(tmp, 'reads')
    for folder in (reads, os.path.join(tmp, 'db'), os.path.join(tmp, 'hostdb'), os.path.join(tmp, 'adapters')):
        os.makedirs(folder, exist_ok=True)
    for sample in samples:
        for mate in (('_1', '_2') if paired else ('',)):
            open(os.path.join(reads, f'{sample}{mate}.fastq.gz'), 'w').close()
    sparse_file(os.path.join(tmp, 'db', 'hash.k2d'), 3000)
    sparse_file(os.path.join(tmp, 'hostdb', 'hash.k2d'), 1000)
    open(os.path.join(tmp, 'db', 'database150mers.kmer_distrib'), 'w').close()
    return {'input_data': {'raw_read_dir': reads, 'paired': paired, 'ext': 'fastq.gz'},
            'krakenDB_dir': os.path.join(tmp, 'db'), 'krakenDB_host': os.path.join(tmp, 'hostdb'),
            'output_dir': os.path.join(tmp, 'output'), 'trimmomatic': {'adapter_dir': os.path.join(tmp, 'adapters')}}


def snakemake(snakefile, config, workdir, *extra, path=None):
    """Run Snakemake on ``config``; returns the completed process (stdout and stderr together)."""
    os.makedirs(workdir, exist_ok=True)
    config_file = os.path.join(workdir, 'config.yaml')
    with open(config_file, 'w') as handle:
        yaml.safe_dump(config, handle)
    environment = dict(os.environ, PATH=path) if path else None
    return subprocess.run([sys.executable, '-m', 'snakemake', '-s', snakefile, '--configfile', config_file,
                           '--directory', workdir, '-c1', *extra],
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=environment, timeout=300)


FAKE_SLURM = r'''#!{python}
import os, sys
name, a, state = os.path.basename(sys.argv[0]), sys.argv[1:], os.environ['FAKE_SLURM_DIR']
def read(file_name):
    path = os.path.join(state, file_name)
    return open(path).read() if os.path.exists(path) else ''
if name in ('sbatch', 'sj-submit'):
    with open(os.path.join(state, 'submitted.txt'), 'w') as f:
        f.write('\n'.join(a) + '\n--- environment\n' + '\n'.join(sorted(os.environ)))
    if read('refuse'):
        sys.stderr.write('sbatch: error: invalid partition specified\n'); sys.exit(1)
    if name == 'sbatch':
        print('Submitted batch job 123')
    else:
        print('Submitted job 456'); sys.stderr.write('Partial success: job 456 is submitted; notifier failed\n'); sys.exit(1)
elif name == 'squeue':
    sys.stdout.write(read('squeue'))
elif name == 'sacct':
    sys.stdout.write(read('sacct'))
elif name == 'scancel':
    open(os.path.join(state, 'cancelled.txt'), 'w').write(' '.join(a))
'''


def fake_slurm(folder, state_dir):
    """Write fake sbatch, sj-submit, squeue, sacct and scancel in ``folder``; they read and write ``state_dir``."""
    os.makedirs(folder, exist_ok=True)
    os.makedirs(state_dir, exist_ok=True)
    for tool in ('sbatch', 'sj-submit', 'squeue', 'sacct', 'scancel'):
        path = os.path.join(folder, tool)
        with open(path, 'w') as handle:
            handle.write(FAKE_SLURM.replace('{python}', sys.executable))
        os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR)


def run_and_measure(module, *arguments):
    """Run ``python -m module arguments`` in a new process; returns (exit code, output, its own peak memory in KB)."""
    code = ('import importlib, resource, sys\n'
            'status = importlib.import_module(sys.argv[1]).main(sys.argv[2:])\n'
            'print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)\n'
            'sys.exit(status)')
    result = subprocess.run([sys.executable, '-c', code, module, *arguments], capture_output=True, text=True)
    lines = result.stdout.strip().splitlines()
    return result.returncode, result.stdout + result.stderr, int(lines[-1]) if lines and lines[-1].isdigit() else 0
