import os
import subprocess
import sys
import tempfile
from unittest import TestCase, mock

from moonstone.reads import WORKFLOW_DIR, launch
from tests.reads.helpers import fake_slurm


class TestLaunch(TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.run_dir = os.path.join(self.tmp, 'run')
        self.state = os.path.join(self.tmp, 'slurm')
        os.makedirs(self.run_dir)
        fake_slurm(os.path.join(self.tmp, 'bin'), self.state)
        patcher = mock.patch.dict(os.environ, {
            'PATH': os.path.join(self.tmp, 'bin') + os.pathsep + '/usr/bin' + os.pathsep + '/bin',
            'FAKE_SLURM_DIR': self.state, 'SLURM_MEM_PER_CPU': '4000', 'SLURM_MEM_PER_NODE': '8000'})
        patcher.start()
        self.addCleanup(patcher.stop)
        self.command = launch.snakemake_command(self.run_dir, profile='/a profile')

    def write(self, name, text, folder=None):
        folder = folder or self.state
        os.makedirs(folder, exist_ok=True)
        with open(os.path.join(folder, name), 'w') as handle:
            handle.write(text)

    def submitted(self):
        with open(os.path.join(self.state, 'submitted.txt')) as handle:
            return handle.read()

    def test_snakemake_command(self):
        local = launch.snakemake_command('/run', cores=8)
        self.assertEqual(local[:3], [sys.executable, '-m', 'snakemake'])
        self.assertEqual(local[3:], ['-s', os.path.join(WORKFLOW_DIR, 'Snakefile'), '--configfile', '/run/config.yaml',
                                     '--directory', '/run', '-c', '8'])
        self.assertEqual(launch.snakemake_command('/run', profile='/p', dry_run=True)[-3:], ['--profile', '/p', '-n'])

    def test_submit_with_sbatch(self):
        job_id = launch.submit(self.run_dir, self.command, head_options='-p atm -c 1')
        self.assertEqual(job_id, '123')
        with open(os.path.join(self.run_dir, 'head_job_id.txt')) as handle:
            self.assertEqual(handle.read().strip(), '123')
        arguments = self.submitted().split('\n--- environment')[0].splitlines()
        self.assertEqual(arguments[:5], ['-p', 'atm', '-c', '1', '--wrap'])
        self.assertIn("--profile '/a profile'", arguments[5])       # the wrapped command is shell-quoted

    def test_submit_with_helper_partial_success(self):
        """The helper prints the job ID, then fails on something else: the job is submitted and recorded."""
        job_id = launch.submit(self.run_dir, self.command, submit='sj-submit --project reads --')
        self.assertEqual(job_id, '456')
        with open(os.path.join(self.run_dir, 'head_job_id.txt')) as handle:
            self.assertEqual(handle.read().strip(), '456')
        self.assertEqual(self.submitted().splitlines()[:3], ['--project', 'reads', '--'])

    def test_submit_failure_reports_the_helper_output(self):
        self.write('refuse', 'yes')
        with self.assertRaisesRegex(RuntimeError, 'invalid partition specified'):
            launch.submit(self.run_dir, self.command)
        self.assertFalse(os.path.exists(os.path.join(self.run_dir, 'head_job_id.txt')))

    def test_submit_drops_inherited_memory_variables(self):
        launch.submit(self.run_dir, self.command)
        environment = self.submitted().split('--- environment\n')[1].splitlines()
        self.assertNotIn('SLURM_MEM_PER_CPU', environment)
        self.assertNotIn('SLURM_MEM_PER_NODE', environment)
        self.assertIn('FAKE_SLURM_DIR', environment)

    def test_status_states(self):
        self.assertEqual(launch.status(self.run_dir).state, 'not started')
        self.write('head_job_id.txt', '123\n', self.run_dir)
        for squeue, sacct, expected in [('PENDING\n', '', 'pending'), ('RUNNING\n', '', 'running'),
                                        ('', 'FAILED\n', 'failed'), ('', 'TIMEOUT\n', 'failed'),
                                        ('', 'CANCELLED by 1000\n', 'cancelled'), ('', 'COMPLETED\n', 'failed'),
                                        ('', '', 'unknown')]:
            with self.subTest(squeue=squeue, sacct=sacct):
                self.write('squeue', squeue)
                self.write('sacct', sacct)
                status = launch.status(self.run_dir)
                self.assertEqual((status.state, status.job_id), (expected, '123'))
        self.write('pipeline_report.md', '# report\n', os.path.join(self.run_dir, 'output'))
        self.assertEqual(launch.status(self.run_dir).state, 'finished')

    def test_status_without_slurm_tools(self):
        self.write('head_job_id.txt', '123\n', self.run_dir)
        with mock.patch.dict(os.environ, {'PATH': '/nonexistent'}):
            self.assertEqual(launch.status(self.run_dir).state, 'unknown')

    def test_status_progress_from_the_snakemake_log(self):
        logs = os.path.join(self.run_dir, '.snakemake', 'log')
        self.write('2026-01-01T000000.snakemake.log', 'old\n3 of 40 steps (8%) done\n', logs)
        self.write('2026-01-02T000000.snakemake.log', 'Building DAG\n5 of 40 steps (12%) done\nrule bracken\n'
                                                      '12 of 40 steps (30%) done\nlast line\n', logs)
        status = launch.status(self.run_dir)
        self.assertEqual((status.done, status.total), (12, 40))
        self.assertTrue(status.log_tail.endswith('last line\n'))

    def test_status_does_not_need_the_submitting_process(self):
        launch.submit(self.run_dir, self.command)
        self.write('squeue', 'RUNNING\n')
        code = f"from moonstone.reads import launch; print(launch.status({self.run_dir!r}).state)"
        result = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True)
        self.assertEqual(result.stdout.strip(), 'running')

    def test_cancel(self):
        with self.assertRaises(ValueError):
            launch.cancel(self.run_dir)
        launch.submit(self.run_dir, self.command)
        launch.cancel(self.run_dir)
        with open(os.path.join(self.state, 'cancelled.txt')) as handle:
            self.assertEqual(handle.read(), '123')
