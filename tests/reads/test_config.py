import os
import tempfile
from unittest import TestCase, mock

import yaml

from moonstone.reads import WORKFLOW_DIR, config
from tests.reads.helpers import fake_bin, make_inputs


class TestConfig(TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.path = fake_bin(os.path.join(self.tmp, 'bin'))
        patcher = mock.patch.dict(os.environ, {'PATH': self.path})
        patcher.start()
        self.addCleanup(patcher.stop)

    def messages(self, settings):
        return [finding.message for finding in config.validate(settings) if finding.level == 'error']

    def test_a_good_configuration_has_no_finding(self):
        self.assertEqual(config.validate(make_inputs(self.tmp)), [])

    def test_build_matches_the_workflow_keys(self):
        built = config.build('/reads', '/out', True, 'fq.gz', {
            'db': '/db', 'host_db': '/hostdb', 'remove_host': ['human', 'mouse'], 'read_len': 100,
            'bracken_threshold': 7, 'keep_trimmed_reads': True, 'min_reads': 5, 'min_classified_pct': 1.5})
        with open(os.path.join(WORKFLOW_DIR, 'config_all_options.yaml')) as handle:
            documented = yaml.safe_load(handle)
        self.assertLessEqual(set(built), set(documented))
        self.assertEqual(built['input_data'], {'raw_read_dir': '/reads', 'paired': True, 'ext': 'fq.gz'})
        self.assertEqual(built['to_perform'], {'trimmomatic': True, 'host_mapping': True, 'bracken': True})
        self.assertEqual((built['krakenDB_dir'], built['krakenDB_host']), ('/db', '/hostdb'))
        self.assertEqual(built['host_taxid'], [9606, 10090])
        self.assertEqual((built['keep_trimmed_reads'], built['keep_kraken_output']), (True, False))
        self.assertEqual(built['bracken'], {'read_len': 100, 'threshold': 7})
        self.assertEqual(built['report'], {'min_reads': 5, 'min_classified_pct': 1.5})
        self.assertEqual(yaml.safe_load(yaml.safe_dump(built)), built)

    def test_build_defaults_equal_moonbeam(self):
        built = config.build('/reads', '/out', True, 'fastq.gz', {'db': '/db', 'host_db': '/hostdb'})
        self.assertEqual(built['host_taxid'], [9606, 10090])
        self.assertEqual(built['bracken'], {'read_len': 150, 'threshold': 5})
        self.assertEqual(built['report'], {'min_reads': 100000, 'min_classified_pct': 20})
        self.assertEqual((built['keep_trimmed_reads'], built['keep_kraken_output']), (False, False))

    def test_hosts_become_taxids(self):
        one = config.build('/r', '/o', True, 'fastq.gz', {'db': '/db', 'host_db': '/h', 'remove_host': ['mouse']})
        self.assertEqual((one['host_taxid'], one['to_perform']['host_mapping']), ([10090], True))
        none = config.build('/r', '/o', True, 'fastq.gz', {'db': '/db', 'host_db': '/h', 'remove_host': []})
        self.assertFalse(none['to_perform']['host_mapping'])
        self.assertNotIn('krakenDB_host', none)
        with self.assertRaisesRegex(ValueError, 'dog'):
            config.build('/r', '/o', True, 'fastq.gz', {'db': '/db', 'remove_host': ['dog']})

    def test_no_host_database_needed_without_host_removal(self):
        settings = make_inputs(self.tmp)
        del settings['krakenDB_host']
        settings['to_perform'] = {'host_mapping': False}
        self.assertEqual(config.validate(settings), [])

    def test_validate_messages(self):
        good = make_inputs(self.tmp)
        cases = [
            ({'input_data': dict(good['input_data'], raw_read_dir='/nonexistent')}, 'not a directory: /nonexistent'),
            ({'input_data': dict(good['input_data'], ext='bam')}, 'No input files found in'),
            ({'krakenDB_dir': os.path.join(self.tmp, 'adapters')}, 'krakenDB_dir is not a Kraken2 database (no hash.k2d)'),
            ({'krakenDB_host': '/nonexistent'}, 'krakenDB_host is not a directory: /nonexistent'),
            ({'bracken': {'read_len': 100}}, 'database100mers.kmer_distrib'),
            ({'krakenDB_dir': ''}, 'Missing required setting: krakenDB_dir'),
            ({'krakenDB_host': ''}, 'Missing required setting: krakenDB_host'),
            ({'keep_intermediates': True}, 'use keep_trimmed_reads and keep_kraken_output'),
        ]
        for change, expected in cases:
            with self.subTest(expected=expected):
                self.assertTrue(any(expected in message for message in self.messages({**good, **change})),
                                self.messages({**good, **change}))

    def test_wrong_mate_names(self):
        settings = make_inputs(self.tmp)
        reads = settings['input_data']['raw_read_dir']
        for name in os.listdir(reads):
            os.rename(os.path.join(reads, name), os.path.join(reads, name.replace('_1.', '_R1.').replace('_2.', '_R2.')))
        self.assertTrue(any('expects {sample}_1/{sample}_2 files' in message for message in self.messages(settings)))

    def test_host_removal_program_must_be_installed(self):
        path = fake_bin(os.path.join(self.tmp, 'bin2'), without=('extract_kraken_reads.py',))
        with mock.patch.dict(os.environ, {'PATH': path}):
            self.assertTrue(any('KrakenTools' in message for message in self.messages(make_inputs(self.tmp))))

    def test_site_file_from_environment_variable(self):
        site_file = os.path.join(self.tmp, 'site.yaml')
        with open(site_file, 'w') as handle:
            handle.write('db: /shared/db\nsubmit: sj-submit --project p --\n')
        with mock.patch.dict(os.environ, {config.SITE_FILE_VARIABLE: site_file}):
            self.assertEqual(config.site(), {'db': '/shared/db', 'submit': 'sj-submit --project p --'})

    def test_no_site_file(self):
        with mock.patch.dict(os.environ, {config.SITE_FILE_VARIABLE: os.path.join(self.tmp, 'absent.yaml')}):
            self.assertEqual(config.site(), {})

    def test_write_never_overwrites(self):
        target = os.path.join(self.tmp, 'config.yaml')
        config.write({'krakenDB_dir': '/db'}, target)
        with open(target) as handle:
            self.assertEqual(yaml.safe_load(handle), {'krakenDB_dir': '/db'})
        with self.assertRaises(FileExistsError):
            config.write({}, target)
