"""Moonstone interface. Start it with the ``moonstone-ui`` command (or: streamlit run moonstone/ui/app.py)."""
import streamlit as st

from moonstone import __version__
from moonstone.ui import views
from moonstone.ui.catalog import TOOLS


def _tool_page(tool):
    def page():
        views.render_tool(tool)
    return st.Page(page, title=tool.title, icon=tool.icon, url_path=tool.key)


st.set_page_config(page_title='Moonstone', page_icon=':material/diamond:', layout='wide')
views.keep_form_values()

pages = {tool.key: _tool_page(tool) for tool in TOOLS.values()}
home = st.Page(lambda: views.render_home(pages), title='Home', icon=':material/home:', url_path='home', default=True)
pages['tutorials'] = st.Page(views.render_tutorials, title='Tutorials', icon=':material/school:',
                             url_path='tutorials')
pages['runs'] = st.Page(views.render_runs, title='Previous runs', icon=':material/history:', url_path='runs')
pages['check'] = st.Page(views.render_check, title='Check my files', icon=':material/fact_check:', url_path='check')
pages['api'] = st.Page(views.render_api, title='Python API', icon=':material/menu_book:', url_path='api')
sections = {
    'Start': [home, pages['tutorials'], pages['runs']],
    'Prepare': [pages['check']] + [pages[tool.key] for tool in TOOLS.values() if tool.group == 'Prepare'],
    'Analyse': [pages[tool.key] for tool in TOOLS.values() if tool.group == 'Analyse'],
    'Reference': [pages['api']],
}
page = st.navigation(sections)
with st.sidebar:
    views.runs_root()
    st.text_input('Run directory', key='runs_root',
                  help='Each run is written to a new directory under this one. Nothing is ever overwritten.')
    st.caption(f'Moonstone {__version__}')
page.run()
