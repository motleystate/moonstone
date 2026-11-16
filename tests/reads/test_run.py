import contextlib
import io
import os
import tempfile
from unittest import TestCase, mock, skipUnless

import yaml

from moonstone.reads import cli, config, run
from tests.reads.helpers import HAS_SNAKEMAKE, fake_bin, fake_slurm, make_inputs


class TestRun(TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.inputs = make_inputs(self.tmp)
        self.reads = self.inputs['input_data']['raw_read_dir']
        self.out = os.path.join(self.tmp, 'run')
        self.state = os.path.join(self.tmp, 'slurm')
        fake_slurm(os.path.join(self.tmp, 'bin'), self.state)
        site_file = os.path.join(self.tmp, 'site.yaml')
        with open(site_file, 'w') as handle:
            yaml.safe_dump({'adapter_dir': os.path.join(self.tmp, 'adapters')}, handle)
        patcher = mock.patch.dict(os.environ, {'PATH': fake_bin(os.path.join(self.tmp, 'bin')),
                                               'FAKE_SLURM_DIR': self.state, config.SITE_FILE_VARIABLE: site_file})
        patcher.start()
        self.addCleanup(patcher.stop)
        self.databases = ['--db', self.inputs['krakenDB_dir'], '--host-db', self.inputs['krakenDB_host']]

    def call(self, *arguments):
        output = io.StringIO()
        with contextlib.redirect_stdout(output), contextlib.redirect_stderr(output):
            try:
                code = cli.main(list(arguments))
            except SystemExit as stop:
                code = stop.code
        return code, output.getvalue()

    def written_config(self):
        with open(os.path.join(self.out, 'config.yaml')) as handle:
            return yaml.safe_load(handle)

    def test_parser_defaults(self):
        args = run.build_parser().parse_args(['reads', 'out', '--db', '/db'])
        self.assertEqual((args.read_len, args.bracken_threshold, args.min_reads, args.min_classified_pct),
                         (150, 5, 100000, 20.0))
        self.assertEqual(args.remove_host, ['human', 'mouse'])
        self.assertFalse(args.keep_trimmed_reads or args.keep_kraken_output or args.slurm or args.single_end)

    def test_help_of_the_keep_switches(self):
        actions = {action.dest: action.help for action in run.build_parser()._actions}
        self.assertIn('allows mapping with another database without trimming again', actions['keep_trimmed_reads'])
        self.assertIn('per-sample reports are always kept', actions['keep_kraken_output'])
        self.assertIn('human (9606)', actions['remove_host'])
        self.assertIn('mouse (10090)', actions['remove_host'])

    def test_site_file_gives_the_database_defaults(self):
        with open(os.environ[config.SITE_FILE_VARIABLE], 'w') as handle:
            handle.write('db: /shared/db\nhost_db: /shared/hosts\n')
        args = run.build_parser().parse_args(['reads', 'out'])
        self.assertEqual((args.db, args.host_db), ('/shared/db', '/shared/hosts'))

    @skipUnless(HAS_SNAKEMAKE, 'needs Snakemake')
    def test_local_dry_run_end_to_end(self):
        code, output = self.call('run', self.reads, self.out, *self.databases, '--dry-run')
        self.assertEqual(code, 0, output)
        self.assertEqual(sorted(os.listdir(os.path.join(self.out, 'reads'))),
                         ['S1_1.fastq.gz', 'S1_2.fastq.gz', 'S2_1.fastq.gz', 'S2_2.fastq.gz'])
        with open(os.path.join(self.out, 'samples.tsv')) as handle:
            self.assertEqual(len(handle.readlines()), 3)
        written = self.written_config()
        self.assertEqual(written['input_data']['raw_read_dir'], os.path.join(self.out, 'reads'))
        self.assertEqual(written['output_dir'], os.path.join(self.out, 'output'))
        self.assertFalse(os.path.exists(os.path.join(self.out, 'output')))      # listed, not run

    @skipUnless(HAS_SNAKEMAKE, 'needs Snakemake')
    def test_local_run_end_to_end(self):
        code, output = self.call('run', self.reads, self.out, *self.databases, '--cores', '2')
        self.assertEqual(code, 0, output)
        self.assertTrue(os.path.exists(os.path.join(self.out, 'output', 'mapping', 'all_samples.bracken.tsv')))
        code, output = self.call('status', self.out)
        self.assertIn('finished', output)

    def test_host_menu(self):
        self.call('run', self.reads, self.out, '--db', self.inputs['krakenDB_dir'], '--host-db',
                  self.inputs['krakenDB_host'], '--remove-host', 'mouse', '--slurm')
        self.assertEqual(self.written_config()['host_taxid'], [10090])
        code, output = self.call('run', self.reads, self.out + '2', *self.databases, '--remove-host', 'dog')
        self.assertEqual(code, 2)
        self.assertIn('invalid choice', output)

    def test_no_host_removal_needs_no_host_database(self):
        code, output = self.call('run', self.reads, self.out, '--db', self.inputs['krakenDB_dir'], '--remove-host',
                                 '--slurm')
        self.assertEqual(code, 0, output)
        written = self.written_config()
        self.assertFalse(written['to_perform']['host_mapping'])
        self.assertNotIn('krakenDB_host', written)

    def test_keep_switches_reach_the_config(self):
        self.call('run', self.reads, self.out, *self.databases, '--keep-trimmed-reads', '--slurm')
        written = self.written_config()
        self.assertEqual((written['keep_trimmed_reads'], written['keep_kraken_output']), (True, False))

    def test_slurm_submission_returns_at_once_with_the_job_id(self):
        code, output = self.call('run', self.reads, self.out, *self.databases, '--slurm')
        self.assertEqual(code, 0, output)
        self.assertIn('123', output)
        with open(os.path.join(self.out, 'head_job_id.txt')) as handle:
            self.assertEqual(handle.read().strip(), '123')
        with open(os.path.join(self.state, 'submitted.txt')) as handle:
            wrapped = handle.read()
        self.assertIn('--profile', wrapped)
        self.assertIn('-m snakemake', wrapped)

    def test_bad_reads_are_reported_before_anything_is_written(self):
        open(os.path.join(self.reads, 'S3_1.fastq.gz'), 'w').close()       # a mate without its pair
        code, output = self.call('run', self.reads, self.out, *self.databases, '--slurm')
        self.assertEqual(code, 1)
        self.assertIn('S3', output)
        self.assertFalse(os.path.exists(self.out))

    def test_bad_database_is_reported_before_anything_is_written(self):
        code, output = self.call('run', self.reads, self.out, '--db', self.tmp, '--host-db', self.inputs['krakenDB_host'])
        self.assertEqual(code, 1)
        self.assertIn('not a Kraken2 database', output)
        self.assertFalse(os.path.exists(self.out))

    def test_a_directory_that_holds_a_run_is_refused(self):
        self.call('run', self.reads, self.out, *self.databases, '--slurm')
        code, output = self.call('run', self.reads, self.out, *self.databases, '--slurm')
        self.assertEqual(code, 1)
        self.assertIn('already holds a run', output)

    def test_other_naming_conventions_and_sub_folders(self):
        nested = os.path.join(self.tmp, 'nested', 'batch1')
        os.makedirs(nested)
        for name in ('A_R1_001.fq.gz', 'A_R2_001.fq.gz'):
            open(os.path.join(nested, name), 'w').close()
        code, output = self.call('run', os.path.join(self.tmp, 'nested'), self.out, *self.databases, '--slurm')
        self.assertEqual(code, 0, output)
        self.assertEqual(sorted(os.listdir(os.path.join(self.out, 'reads'))), ['A_1.fq.gz', 'A_2.fq.gz'])
        self.assertEqual(self.written_config()['input_data']['ext'], 'fq.gz')

    def test_status_and_cancel_commands(self):
        self.call('run', self.reads, self.out, *self.databases, '--slurm')
        with open(os.path.join(self.state, 'squeue'), 'w') as handle:
            handle.write('RUNNING\n')
        code, output = self.call('status', self.out)
        self.assertEqual(code, 0)
        self.assertIn('running', output)
        self.assertIn('123', output)
        self.assertEqual(self.call('cancel', self.out)[0], 0)
        self.assertTrue(os.path.exists(os.path.join(self.state, 'cancelled.txt')))

    def test_cancel_without_slurm_is_an_error_message_not_a_traceback(self):
        self.call('run', self.reads, self.out, *self.databases, '--slurm')
        with mock.patch.dict(os.environ, {'PATH': '/nonexistent'}):
            code, output = self.call('cancel', self.out)
        self.assertEqual(code, 1)
        self.assertIn('scancel', output)

    def test_usage(self):
        code, output = self.call()
        self.assertIn('run', output)
        self.assertIn('status', output)
        self.assertEqual(self.call('nonsense')[0], 2)


class TestCheckInputs(TestCase):
    """run.check_inputs is what the interface calls: the checks of the command, without writing anything."""

    def test_same_findings_as_the_command_and_nothing_written(self):
        tmp = tempfile.mkdtemp()
        inputs = make_inputs(tmp)
        reads = inputs['input_data']['raw_read_dir']
        with mock.patch.dict(os.environ, {'PATH': fake_bin(os.path.join(tmp, 'bin')),
                                          config.SITE_FILE_VARIABLE: os.path.join(tmp, 'absent.yaml')}):
            options = {'db': inputs['krakenDB_dir'], 'host_db': inputs['krakenDB_host'], 'adapter_dir': None}
            found, findings = run.check_inputs(reads, options)
            self.assertEqual([sample.name for sample in found], ['S1', 'S2'])
            self.assertTrue(all('adapter' in finding.message for finding in findings), findings)
            found, findings = run.check_inputs(reads, dict(options, db=None))
            self.assertIn('--db', findings[0].message)
        self.assertEqual(sorted(os.listdir(reads)), ['S1_1.fastq.gz', 'S1_2.fastq.gz', 'S2_1.fastq.gz', 'S2_2.fastq.gz'])
