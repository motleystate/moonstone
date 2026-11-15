"""The ``moonstone-ui`` command: start the interface in the browser."""
import os
import subprocess
import sys


def main(argv=None) -> int:
    """Run the Streamlit app on this machine only; extra arguments are passed to ``streamlit run``."""
    try:
        import streamlit  # noqa: F401
    except ImportError:
        print('The interface needs Streamlit: pip install "moonstone[ui]"')
        return 1
    app = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'app.py')
    command = [sys.executable, '-m', 'streamlit', 'run', app, '--server.address', '127.0.0.1',
               '--browser.gatherUsageStats', 'false',
               # one accent colour per theme, so the app still follows the system's light or dark setting
               '--theme.light.primaryColor', '#4F6FA8', '--theme.dark.primaryColor', '#8FB0E6']
    return subprocess.call(command + list(sys.argv[1:] if argv is None else argv))


if __name__ == '__main__':
    sys.exit(main())
