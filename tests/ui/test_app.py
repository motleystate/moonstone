import os
import subprocess
import sys
import tempfile
from unittest import TestCase

import pandas as pd
from streamlit.testing.v1 import AppTest

from moonstone import examples
from moonstone.main import build_parser
import moonstone.ui
from moonstone.ui.catalog import TOOLS

TIMEOUT = 120
APP = os.path.join(os.path.dirname(moonstone.ui.__file__), 'app.py')    # not imported: importing runs it


def tool_page(key, runs_root=None):
    import streamlit as st
    from moonstone.ui import catalog, views
    if runs_root:
        st.session_state.setdefault('runs_root', runs_root)
    views.render_tool(catalog.TOOLS[key])


def results_page(run_dir):
    from pathlib import Path
    from moonstone.ui import views
    views.render_results(Path(run_dir))


class TestApp(TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def page(self, key):
        return AppTest.from_function(tool_page, args=(key, self.tmp), default_timeout=TIMEOUT).run()

    def test_app_starts_on_a_home_page_with_one_card_per_tool(self):
        at = AppTest.from_file(APP, default_timeout=TIMEOUT).run()
        self.assertFalse(at.exception)
        for key in TOOLS:
            self.assertIsNotNone(at.button(key=f'home.{key}'))

    def test_every_tool_page_renders_empty(self):
        for key in TOOLS:
            with self.subTest(tool=key):
                at = self.page(key)
                self.assertFalse(at.exception)
                self.assertTrue(at.button(key=f'{key}.run').disabled)

    def test_example_button_fills_a_valid_form(self):
        at = self.page('permanova')
        at.button(key='permanova.example').click().run()
        self.assertFalse(at.exception)
        self.assertEqual(at.multiselect(key='permanova.permanova').value, ['disease'])
        self.assertEqual(at.selectbox(key='permanova.strata').value, 'batch')
        self.assertFalse(at.button(key='permanova.run').disabled)
        self.assertFalse(at.error)
        self.assertIn('--permanova disease --covariates age antibiotics --strata batch', at.code[0].value)

    def test_every_example_makes_a_valid_form(self):
        for key in TOOLS:
            if not TOOLS[key].example:
                continue
            with self.subTest(tool=key):
                at = self.page(key)
                at.button(key=f'{key}.example').click().run()
                self.assertFalse(at.exception)
                self.assertFalse(at.error)
                self.assertFalse(at.button(key=f'{key}.run').disabled)

    def test_help_text_comes_from_the_parser(self):
        at = self.page('permanova')
        at.button(key='permanova.example').click().run()
        actions = {action.dest: action for action in build_parser()._actions}
        self.assertEqual(at.selectbox(key='permanova.strata').help, actions['strata'].help)
        self.assertEqual(at.selectbox(key='permanova.taxa_level').help, actions['taxa_level'].help)
        self.assertEqual(at.selectbox(key='permanova.taxa_level').value, 'species')

    def test_command_error_disables_run(self):
        at = self.page('correlation')
        at.button(key='correlation.example').click().run()
        at.number_input(key='correlation.alpha').set_value(2.0).run()
        self.assertIn('--alpha must be between 0 and 1.', [error.value for error in at.error])
        self.assertTrue(at.button(key='correlation.run').disabled)

    def test_unreadable_metadata_keeps_the_form(self):
        at = self.page('permanova')
        at.text_input(key='permanova.countfile').set_value(examples.path('species_counts.tsv'))
        at.text_input(key='permanova.metadata').set_value('/nope.csv').run()
        self.assertFalse(at.exception)
        self.assertTrue(any('/nope.csv' in error.value for error in at.error))
        self.assertTrue(at.button(key='permanova.run').disabled)
        self.assertIsNotNone(at.text_input(key='permanova.permanova.text'))

    def test_wrong_table_is_reported_with_the_reader_message(self):
        at = self.page('permanova')
        at.text_input(key='permanova.countfile').set_value(examples.path('resistance_genes.csv'))
        at.text_input(key='permanova.metadata').set_value(examples.path('metadata.csv')).run()
        self.assertTrue(any('no taxonomy column found' in error.value for error in at.error))

    def test_run_shows_results(self):
        at = self.page('read_threshold')
        at.button(key='read_threshold.example').click().run()
        at.button(key='read_threshold.run').click().run()
        self.assertFalse(at.exception)
        self.assertFalse(at.error)
        self.assertTrue(len(at.dataframe) >= 1)
        self.assertEqual(len(os.listdir(self.tmp)), 1)

    def test_failed_run_shows_the_command_line_message(self):
        at = self.page('permanova')
        at.button(key='permanova.example').click().run()
        # typed in the text fallback of an unknown metadata: a column that does not exist
        at.selectbox(key='permanova.strata').set_value(None)
        at.multiselect(key='permanova.covariates').set_value([])
        at.multiselect(key='permanova.select').set_value([]).run()
        at.selectbox(key='permanova.group_column').set_value('Type').run()   # two samples per Type: refused
        at.button(key='permanova.run').click().run()
        self.assertFalse(at.exception)
        self.assertTrue(at.error)
        self.assertEqual(len(os.listdir(self.tmp)), 1)

    def test_results_preview_is_capped(self):
        run_dir = os.path.join(self.tmp, '20260101-000000-permanova')
        os.mkdir(run_dir)
        pd.DataFrame({'x': range(5000)}).to_csv(os.path.join(run_dir, 'permanova.csv'), index=False)
        with open(os.path.join(run_dir, 'command.txt'), 'w') as handle:
            handle.write('moonstone ...\n')
        at = AppTest.from_function(results_page, args=(run_dir,), default_timeout=TIMEOUT).run()
        self.assertFalse(at.exception)
        self.assertEqual(len(at.dataframe[0].value), 1000)
        self.assertTrue(any('5,000' in caption.value for caption in at.caption))

    def test_launcher_without_streamlit(self):
        code = "import sys; sys.modules['streamlit'] = None; from moonstone.ui.launch import main; sys.exit(main([]))"
        result = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True)
        self.assertEqual(result.returncode, 1)
        self.assertIn('pip install "moonstone[ui]"', result.stdout)


def check_page():
    from moonstone.ui import views
    views.render_check()


def runs_page(runs_root):
    import streamlit as st
    from moonstone.ui import views
    st.session_state.setdefault('runs_root', runs_root)
    views.render_runs()


def api_page():
    from moonstone.ui import views
    views.render_api()


class TestOtherPages(TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def test_check_page_offers_a_conversion_and_never_overwrites(self):
        source = os.path.join(self.tmp, 'genes.csv')
        with open(source, 'w') as handle:
            handle.write('gene;S001;S002\ntetM;1,5;2,25\n')
        at = AppTest.from_function(check_page, default_timeout=TIMEOUT).run()
        at.text_input(key='check.file').set_value(source).run()
        self.assertFalse(at.exception)
        self.assertTrue(any('decimal mark' in error.value for error in at.error))
        at.button(key='check.convert').click().run()
        self.assertFalse(at.exception)
        converted = os.path.join(self.tmp, 'genes.moonstone.csv')
        with open(converted) as handle:
            self.assertEqual(handle.read(), 'gene,S001,S002\ntetM,1.5,2.25\n')
        with open(source) as handle:
            self.assertEqual(handle.read(), 'gene;S001;S002\ntetM;1,5;2,25\n')
        at.button(key='check.convert').click().run()
        self.assertTrue(any('already exists; it is never overwritten' in error.value for error in at.error))

    def test_clean_file_offers_no_conversion(self):
        at = AppTest.from_function(check_page, default_timeout=TIMEOUT).run()
        at.text_input(key='check.file').set_value(examples.path('metadata.csv')).run()
        self.assertFalse(at.exception)
        self.assertFalse(at.error)
        self.assertTrue(at.success)
        self.assertNotIn('check.convert', [button.key for button in at.button])

    def test_check_page_accepts_a_spreadsheet(self):
        source = os.path.join(self.tmp, 'metadata.xlsx')
        pd.DataFrame({'sample': ['S1', 'S2'], 'age': [25, 31]}).to_excel(source, index=False)
        at = AppTest.from_function(check_page, default_timeout=TIMEOUT).run()
        at.text_input(key='check.file').set_value(source).run()
        self.assertFalse(at.exception)
        self.assertTrue(at.success)

    def test_check_page_for_a_tool_uses_its_reader(self):
        at = AppTest.from_function(check_page, default_timeout=TIMEOUT).run()
        at.text_input(key='check.file').set_value(examples.path('resistance_genes.csv'))
        at.selectbox(key='check.tool').set_value('permanova').run()
        self.assertTrue(any('no taxonomy column found' in error.value for error in at.error))

    def test_previous_runs_lists_finished_and_failed(self):
        for name, code in [('20260101-000000-permanova', '0'), ('20260102-000000-correlation', '1')]:
            os.mkdir(os.path.join(self.tmp, name))
            for file_name, text in [('command.txt', f'moonstone {name}\n'), ('exit_code.txt', code),
                                    ('run.log', 'the log\n')]:
                with open(os.path.join(self.tmp, name, file_name), 'w') as handle:
                    handle.write(text)
        at = AppTest.from_function(runs_page, args=(self.tmp,), default_timeout=TIMEOUT).run()
        self.assertFalse(at.exception)
        self.assertEqual(list(at.dataframe[0].value['Status']), ['failed', 'finished'])
        self.assertTrue(at.error)    # the newest run, opened by default, failed
        at.selectbox(key='runs.pick').set_value('20260101-000000-permanova').run()
        self.assertFalse(at.error)

    def test_no_previous_run(self):
        at = AppTest.from_function(runs_page, args=(self.tmp,), default_timeout=TIMEOUT).run()
        self.assertFalse(at.exception)
        self.assertTrue(at.info)

    def test_api_search(self):
        at = AppTest.from_function(api_page, default_timeout=TIMEOUT).run()
        self.assertFalse(at.exception)
        everything = len(at.expander)
        at.text_input(key='api.search').set_value('Kraken2Parser').run()
        self.assertEqual([expander.label.split('**')[1] for expander in at.expander], ['SunbeamKraken2Parser'])
        self.assertGreater(everything, 60)


def tutorials_page(pick=None):
    import streamlit as st
    from moonstone.ui import views
    if pick:
        st.session_state.setdefault('tutorials.pick', pick)
    views.render_tutorials()


class TestTutorialPages(TestCase):

    def test_tutorials_page_opens_on_the_first_tutorial(self):
        at = AppTest.from_function(tutorials_page, default_timeout=TIMEOUT).run()
        self.assertFalse(at.exception)
        self.assertFalse(at.error)
        self.assertEqual(at.header[0].value, 'Your first analysis, start to finish')
        self.assertGreaterEqual(len(at.dataframe), 3)

    def test_picking_a_tutorial_shows_its_code_and_results(self):
        at = AppTest.from_function(tutorials_page, default_timeout=TIMEOUT).run()
        at.radio(key='tutorials.menu.Describe').set_value('plots').run()
        self.assertFalse(at.exception)
        self.assertEqual(at.header[0].value, 'Plotting counts and metadata')
        self.assertTrue(any('PlotTaxonomyCounts(counts)' in code.value for code in at.code))
        self.assertGreaterEqual(len(at.get('plotly_chart')), 9)

    def test_every_tutorial_renders(self):
        from moonstone.ui.tutorials import TUTORIALS
        for key in TUTORIALS:
            with self.subTest(tutorial=key):
                at = AppTest.from_function(tutorials_page, args=(key,), default_timeout=TIMEOUT).run()
                self.assertFalse(at.exception)
                self.assertFalse(at.error)

    def test_tool_page_has_its_tutorial(self):
        at = AppTest.from_function(tool_page, args=('read_threshold',), default_timeout=TIMEOUT)
        at.session_state['read_threshold.view'] = 'Tutorial'
        at.run()
        self.assertFalse(at.exception)
        self.assertTrue(any('analyze_normalized_reads' in code.value for code in at.code))
        self.assertEqual(len(at.get('image')), 2)


def reads_without_snakemake():
    import sys
    from moonstone.ui import catalog, views
    saved = sys.modules.get('snakemake')
    sys.modules['snakemake'] = None
    try:
        views.render_tool(catalog.TOOLS['reads'])
    finally:
        if saved is None:
            sys.modules.pop('snakemake', None)
        else:
            sys.modules['snakemake'] = saved


class TestReadsPages(TestCase):

    def setUp(self):
        from unittest import mock
        from tests.reads.helpers import fake_bin, fake_slurm, make_inputs
        self.tmp = tempfile.mkdtemp()
        self.runs = os.path.join(self.tmp, 'runs')
        self.inputs = make_inputs(self.tmp)
        self.reads = self.inputs['input_data']['raw_read_dir']
        self.state = os.path.join(self.tmp, 'slurm')
        fake_slurm(os.path.join(self.tmp, 'bin'), self.state)
        with open(os.path.join(self.tmp, 'site.yaml'), 'w') as handle:
            handle.write(f"adapter_dir: {os.path.join(self.tmp, 'adapters')}\n")
        patcher = mock.patch.dict(os.environ, {'PATH': fake_bin(os.path.join(self.tmp, 'bin')),
                                               'FAKE_SLURM_DIR': self.state,
                                               'MOONSTONE_SITE': os.path.join(self.tmp, 'site.yaml')})
        patcher.start()
        self.addCleanup(patcher.stop)

    def page(self, key):
        return AppTest.from_function(tool_page, args=(key, self.runs), default_timeout=TIMEOUT).run()

    def filled_reads_page(self):
        at = self.page('reads')
        at.text_input(key='reads.countfile').set_value(self.reads)
        at.text_input(key='reads.db').set_value(self.inputs['krakenDB_dir'])
        at.text_input(key='reads.host_db').set_value(self.inputs['krakenDB_host']).run()
        return at

    def submitted_run(self, squeue='RUNNING\n', log='Building DAG\n12 of 40 steps (30%) done\n'):
        run_dir = os.path.join(self.runs, '20260101-000000-reads')
        os.makedirs(os.path.join(run_dir, '.snakemake', 'log'))
        for name, text in [('command.txt', 'python -m moonstone.reads.run ...\n'), ('exit_code.txt', '0\n'),
                           ('run.log', 'Submitted to Slurm as job 123.\n'), ('head_job_id.txt', '123\n'),
                           ('.snakemake/log/2026-01-01T000000.snakemake.log', log)]:
            with open(os.path.join(run_dir, name), 'w') as handle:
                handle.write(text)
        with open(os.path.join(self.state, 'squeue'), 'w') as handle:
            handle.write(squeue)
        return run_dir

    def text(self, at):
        return '\n'.join(element.value for element in list(at.markdown) + list(at.caption))

    def test_reads_page_renders_empty_without_an_example_button(self):
        at = self.page('reads')
        self.assertFalse(at.exception)
        self.assertTrue(at.button(key='reads.run').disabled)
        self.assertNotIn('reads.example', [button.key for button in at.button])

    def test_reads_page_with_valid_inputs(self):
        at = self.filled_reads_page()
        self.assertFalse(at.exception)
        self.assertFalse(at.error)
        self.assertFalse(at.button(key='reads.run').disabled)
        self.assertTrue(at.code[0].value.startswith('python -m moonstone.reads.run ' + self.reads))
        self.assertIn('--remove-host human mouse', at.code[0].value)
        self.assertIn('2 samples, paired-end: S1, S2', self.text(at))

    def test_reads_findings_block_the_run(self):
        open(os.path.join(self.reads, 'S3_1.fastq.gz'), 'w').close()
        at = self.filled_reads_page()
        self.assertTrue(any('S3' in error.value for error in at.error))
        self.assertTrue(at.button(key='reads.run').disabled)

    def test_host_menu_shows_names_and_taxids(self):
        at = self.filled_reads_page()
        menu = at.multiselect(key='reads.remove_host')
        self.assertEqual(menu.options, ['Human (9606)', 'Mouse (10090)'])
        self.assertEqual(menu.value, ['human', 'mouse'])
        at.text_input(key='reads.host_db').set_value('')
        menu.set_value([]).run()
        self.assertFalse(at.error)                                   # no host removal: no host database needed
        self.assertRegex(at.code[0].value, r'--remove-host($| --)')    # the flag is sent with no name

    def test_keep_switch_help_comes_from_the_parser(self):
        at = self.page('reads')
        self.assertIn('without trimming again', at.checkbox(key='reads.keep_trimmed_reads').help)
        self.assertIn('per-sample reports are always kept', at.checkbox(key='reads.keep_kraken_output').help)

    def test_slurm_submission_from_the_page(self):
        at = self.filled_reads_page()
        at.checkbox(key='reads.slurm').check().run()
        with open(os.path.join(self.state, 'squeue'), 'w') as handle:
            handle.write('PENDING\n')
        at.button(key='reads.run').click().run()
        self.assertFalse(at.exception)
        self.assertFalse(at.error)
        run_dir = os.path.join(self.runs, os.listdir(self.runs)[0])
        with open(os.path.join(run_dir, 'head_job_id.txt')) as handle:
            self.assertEqual(handle.read().strip(), '123')
        self.assertIn('pending', self.text(at))
        self.assertIn('Slurm job 123', self.text(at))

    def test_submitted_run_shows_status_after_restart(self):
        run_dir = self.submitted_run()
        at = AppTest.from_function(results_page, args=(run_dir,), default_timeout=TIMEOUT).run()
        self.assertFalse(at.exception)
        self.assertIn('running', self.text(at))
        self.assertIn('12 of 40 steps', self.text(at))

    def test_cancel_asks_for_confirmation(self):
        run_dir = self.submitted_run()
        at = AppTest.from_function(results_page, args=(run_dir,), default_timeout=TIMEOUT).run()
        confirm = at.button(key='cancel.20260101-000000-reads')
        self.assertIn('Yes, cancel', confirm.label)
        self.assertFalse(os.path.exists(os.path.join(self.state, 'cancelled.txt')))
        confirm.click().run()
        self.assertFalse(at.exception)
        with open(os.path.join(self.state, 'cancelled.txt')) as handle:
            self.assertEqual(handle.read(), '123')

    def test_a_finished_run_has_no_cancel_button(self):
        run_dir = self.submitted_run(squeue='')
        os.makedirs(os.path.join(run_dir, 'output', 'mapping'))
        with open(os.path.join(run_dir, 'output', 'pipeline_report.md'), 'w') as handle:
            handle.write('# Moonstone reads pipeline report\n\n| Sample | Reads |\n|---|---|\n| S1 | 1000 |\n')
        with open(os.path.join(run_dir, 'output', 'mapping', 'all_samples.bracken.tsv'), 'w') as handle:
            handle.write('# Constructed from biom file\n#OTU ID\tS1.brackenCorrected_taxa\ttaxonomy\n'
                         '562\t600.0\tk__Bacteria\n')
        at = AppTest.from_function(results_page, args=(run_dir,), default_timeout=TIMEOUT).run()
        self.assertFalse(at.exception)
        self.assertIn('finished', self.text(at))
        self.assertNotIn('cancel.20260101-000000-reads', [button.key for button in at.button])
        self.assertIn('Moonstone reads pipeline report', self.text(at))                # the report is rendered
        self.assertEqual(list(at.dataframe[0].value.columns), ['#OTU ID', 'S1.brackenCorrected_taxa', 'taxonomy'])
        at.button(key='next.20260101-000000-reads').click().run()                       # "Make a Moonstone table"
        self.assertEqual(at.session_state['reads_table.countfile'],
                         os.path.join(run_dir, 'output', 'mapping', 'all_samples.bracken.tsv'))

    def test_table_tool_runs_and_offers_normalization(self):
        at = self.page('reads_table')
        at.button(key='reads_table.example').click().run()
        at.button(key='reads_table.run').click().run()
        self.assertFalse(at.exception)
        self.assertFalse(at.error)
        run_dir = os.path.join(self.runs, os.listdir(self.runs)[0])
        output = os.path.join(run_dir, 'counts_raw.csv')
        self.assertTrue(os.path.exists(output))
        at.button(key=f'next.{os.path.basename(run_dir)}').click().run()
        self.assertEqual(at.session_state['normalize.countfile'], output)

    def test_normalize_tool_needs_a_method(self):
        at = self.page('normalize')
        at.text_input(key='normalize.countfile').set_value(examples.path('taxa_counts.csv')).run()
        self.assertTrue(at.button(key='normalize.run').disabled)
        self.assertIn('--method is required.', [error.value for error in at.error])
        at.selectbox(key='normalize.method').set_value('relative_abundance').run()
        self.assertFalse(at.error)
        at.button(key='normalize.run').click().run()
        self.assertFalse(at.exception)
        run_dir = os.path.join(self.runs, os.listdir(self.runs)[0])
        self.assertTrue(os.path.exists(os.path.join(run_dir, 'counts_relative_abundance.csv')))
        self.assertIn('not for PERMANOVA', '\n'.join(code.value for code in at.code))

    def test_wide_table_preview_is_capped(self):
        run_dir = os.path.join(self.runs, '20260101-000000-reads_table')
        os.makedirs(run_dir)
        pd.DataFrame([[1] * 2000], columns=[f'S{i}' for i in range(2000)]).to_csv(
            os.path.join(run_dir, 'counts_raw.csv'), index=False)
        with open(os.path.join(run_dir, 'command.txt'), 'w') as handle:
            handle.write('moonstone-reads table ...\n')
        at = AppTest.from_function(results_page, args=(run_dir,), default_timeout=TIMEOUT).run()
        self.assertFalse(at.exception)
        self.assertEqual(at.dataframe[0].value.shape[1], 50)
        self.assertTrue(any('2,000' in caption.value for caption in at.caption))

    def test_reads_page_without_the_reads_extra(self):
        at = AppTest.from_function(reads_without_snakemake, default_timeout=TIMEOUT).run()
        self.assertFalse(at.exception)
        self.assertTrue(any('moonstone[reads]' in warning.value for warning in at.warning))
        self.assertTrue(at.button(key='reads.run').disabled)

    def test_previous_runs_show_the_slurm_state(self):
        self.submitted_run()
        at = AppTest.from_function(runs_page, args=(self.runs,), default_timeout=TIMEOUT).run()
        self.assertFalse(at.exception)
        self.assertEqual(list(at.dataframe[0].value['Status']), ['running'])
