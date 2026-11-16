import glob
import os
import re
import tempfile
from unittest import TestCase, skipUnless

from tests.reads.helpers import HAS_SNAKEMAKE, fake_bin, make_inputs, snakemake


def job_counts(log):
    """{rule: number of jobs} from the 'Job stats' table of a Snakemake log."""
    table = log.split('Job stats:')[1].split('\n\n')[0]
    return {rule: int(count) for rule, count in re.findall(r'^(\w+)\s+(\d+)\s*$', table, flags=re.M) if rule != 'total'}


@skipUnless(HAS_SNAKEMAKE, 'needs Snakemake (pip install "moonstone[reads]")')
class TestWorkflow(TestCase):

    def setUp(self):
        from moonstone.reads import WORKFLOW_DIR
        self.snakefile = os.path.join(WORKFLOW_DIR, 'Snakefile')
        self.tmp = tempfile.mkdtemp()
        self.path = fake_bin(os.path.join(self.tmp, 'bin'))

    def dry_run(self, config, path=None):
        return snakemake(self.snakefile, config, os.path.join(self.tmp, 'work'), '-n', path=path or self.path)

    def run_workflow(self, config):
        result = snakemake(self.snakefile, config, os.path.join(self.tmp, 'work'), path=self.path)
        self.assertEqual(result.returncode, 0, result.stdout)
        return config['output_dir']

    def test_dry_run_paired_and_single(self):
        for paired in (True, False):
            with self.subTest(paired=paired):
                result = self.dry_run(make_inputs(os.path.join(self.tmp, str(paired)), paired=paired))
                self.assertEqual(result.returncode, 0, result.stdout)
                self.assertIn('trimmomatic', result.stdout)
                # Kraken2 memory = database size + 4000 MB
                self.assertIn('mem_mb=7000', result.stdout)
                self.assertIn('mem_mb=5000', result.stdout)

    def test_missing_read_directory_is_rejected(self):
        config = make_inputs(self.tmp)
        config['input_data']['raw_read_dir'] = '/nonexistent'
        self.assertIn('not a directory', self.dry_run(config).stdout)

    def test_non_database_is_rejected(self):
        config = make_inputs(self.tmp)
        config['krakenDB_dir'] = os.path.join(self.tmp, 'adapters')
        self.assertIn('not a Kraken2 database', self.dry_run(config).stdout)

    def test_missing_bracken_kmer_file_is_rejected(self):
        config = dict(make_inputs(self.tmp), bracken={'read_len': 100})
        self.assertIn('database100mers.kmer_distrib', self.dry_run(config).stdout)

    def test_missing_required_setting_is_named(self):
        config = dict(make_inputs(self.tmp), krakenDB_dir='')
        self.assertIn('Missing required setting: krakenDB_dir', self.dry_run(config).stdout)

    def test_no_host_database_needed_without_host_mapping(self):
        config = dict(make_inputs(self.tmp), krakenDB_host='', to_perform={'host_mapping': False})
        self.assertEqual(self.dry_run(config).returncode, 0)

    def test_old_keep_intermediates_setting_is_rejected(self):
        result = self.dry_run(dict(make_inputs(self.tmp), keep_intermediates=True))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('keep_trimmed_reads', result.stdout)
        self.assertIn('keep_kraken_output', result.stdout)

    def test_host_removal_script_must_be_on_path(self):
        path = fake_bin(os.path.join(self.tmp, 'bin2'), without=('extract_kraken_reads.py',))
        result = self.dry_run(make_inputs(self.tmp), path=path)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('KrakenTools', result.stdout)

    def test_missing_adapters_are_reported(self):
        config = make_inputs(self.tmp)
        del config['trimmomatic']
        result = self.dry_run(config)
        if result.returncode:    # passes where the conda trimmomatic package provides the adapters
            self.assertIn('adapter_dir', result.stdout)

    def test_whole_workflow_with_default_switches(self):
        """Large intermediates are gone; everything a later step or the user needs is there."""
        out = self.run_workflow(dict(make_inputs(self.tmp), report={'min_reads': 100}))
        for kept in ['pipeline_report.md', 'mapping/all_samples.bracken.tsv', 'mapping/all_samples.kraken.tsv',
                     'mapping/all_samples.bracken.biom', 'mapping/kraken2_final/S1.taxa.report',
                     'mapping/kraken2_host/S1.taxa.report', 'mapping/bracken_results/S2.brackenCorrected_taxa.report']:
            self.assertTrue(os.path.exists(os.path.join(out, kept)), kept)
        self.assertEqual(glob.glob(out + '/**/*.mapping.kraken', recursive=True), [])     # per-read Kraken2 files
        self.assertEqual(glob.glob(out + '/**/*.fastq*', recursive=True), [])             # trimmed and host-free reads
        report = open(os.path.join(out, 'pipeline_report.md')).read()
        self.assertIn('| S1 | 1000 | 30.0 | 1000 | 90.0 | 10.0 | 2 |', report)
        self.assertIn('None.', report.split('## Warnings')[1].split('## Output')[0])
        self.assertIn('keep_trimmed_reads', report)

    def test_keep_trimmed_reads(self):
        out = self.run_workflow(dict(make_inputs(self.tmp), keep_trimmed_reads=True))
        kept = sorted(os.path.basename(p) for p in glob.glob(out + '/**/*.fastq*', recursive=True))
        self.assertEqual(kept, ['S1_1.trimmed.fastq.gz', 'S1_2.trimmed.fastq.gz', 'S2_1.trimmed.fastq.gz',
                                'S2_2.trimmed.fastq.gz'])
        self.assertEqual(glob.glob(out + '/**/*.mapping.kraken', recursive=True), [])

    def test_keep_kraken_output(self):
        out = self.run_workflow(dict(make_inputs(self.tmp), keep_kraken_output=True))
        kept = sorted(os.path.relpath(p, out) for p in glob.glob(out + '/**/*.mapping.kraken', recursive=True))
        self.assertEqual(kept, ['mapping/kraken2_final/S1.mapping.kraken', 'mapping/kraken2_final/S2.mapping.kraken',
                                'mapping/kraken2_host/S1.mapping.kraken', 'mapping/kraken2_host/S2.mapping.kraken'])
        self.assertEqual(glob.glob(out + '/**/*.fastq*', recursive=True), [])

    def test_single_end_without_host_removal_or_bracken(self):
        config = dict(make_inputs(self.tmp, paired=False), krakenDB_host='',
                      to_perform={'host_mapping': False, 'bracken': False})
        out = self.run_workflow(config)
        self.assertTrue(os.path.exists(os.path.join(out, 'mapping/all_samples.kraken.tsv')))
        self.assertFalse(os.path.exists(os.path.join(out, 'mapping/all_samples.bracken.tsv')))

    @skipUnless(os.environ.get('MOONBEAM_REPO'), 'set MOONBEAM_REPO to a Moonbeam checkout to compare with it')
    def test_same_jobs_as_moonbeam(self):
        ours = dict(make_inputs(self.tmp), keep_trimmed_reads=True, keep_kraken_output=True)
        theirs = dict(make_inputs(self.tmp), keep_intermediates=True)
        moonbeam = snakemake(os.path.join(os.environ['MOONBEAM_REPO'], 'Snakefile'), theirs,
                             os.path.join(self.tmp, 'moonbeam'), '-n', path=self.path)
        self.assertEqual(moonbeam.returncode, 0, moonbeam.stdout)
        self.assertEqual(job_counts(self.dry_run(ours).stdout), job_counts(moonbeam.stdout))
