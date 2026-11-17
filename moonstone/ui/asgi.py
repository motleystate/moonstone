"""
The interface as an ASGI application, for ``moonstone-ui --remote``.

Streamlit's own server cannot listen on a Unix socket; uvicorn can, and serves this application:
    python -m uvicorn moonstone.ui.asgi:app --uds /path/to/private/folder/ui.sock
"""
import os

import streamlit as st
from streamlit import config

# the settings that `moonstone-ui` passes to `streamlit run` on the command line
for _option, _value in {'browser.gatherUsageStats': False, 'theme.light.primaryColor': '#4F6FA8',
                        'theme.dark.primaryColor': '#8FB0E6'}.items():
    config.set_option(_option, _value)

app = st.App(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'app.py'))
