"""
The ``moonstone-ui`` command: start the interface.

    moonstone-ui             on this computer, in your browser
    moonstone-ui --remote    on a server or cluster: no network port is opened; the interface listens on a Unix
                             socket that only you can open, and you reach it through an SSH tunnel
    moonstone-ui --where     how to connect to the interface started with --remote
"""
import argparse
import os
import signal
import socket
import stat
import subprocess
import sys
import time

PRIVATE_FOLDER = os.path.join('~', '.moonstone')


def default_socket() -> str:
    """~/.moonstone/ui-<host>.sock: one per machine, so servers on two nodes sharing a home do not collide."""
    return os.path.join(os.path.expanduser(PRIVATE_FOLDER), f'ui-{socket.gethostname()}.sock')


def connection_file() -> str:
    return os.path.join(os.path.expanduser(PRIVATE_FOLDER), 'connection.txt')


def prepare_socket(path: str) -> None:
    """
    Make sure the socket can be created where only its owner can reach it.

    Other users of a shared machine can connect to any network port, even one bound to 127.0.0.1, but not to a
    socket inside a folder they cannot enter. The folder is created with mode 0700, and refused if it is open
    to anyone else. A socket left by a server that has stopped is removed.

    :raise RuntimeError: the folder is not private, or a server is already listening on the socket.
    """
    folder = os.path.dirname(os.path.abspath(path))
    os.makedirs(folder, mode=0o700, exist_ok=True)
    status = os.stat(folder)
    if status.st_uid != os.getuid() or stat.S_IMODE(status.st_mode) & 0o077:
        raise RuntimeError(f'{folder} can be entered by other users (mode {stat.S_IMODE(status.st_mode):04o}). '
                           f'Make it private first: chmod 700 {folder}')
    if os.path.exists(path):
        probe = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        try:
            probe.connect(path)
        except OSError:
            os.remove(path)          # left by a server that has stopped
        else:
            raise RuntimeError(f'Moonstone is already serving on {path}. See: moonstone-ui --where')
        finally:
            probe.close()


def connection_text(socket_path: str, port: int = 8501, environ=None, host: str = None) -> str:
    """What to run on your own computer to reach the interface."""
    environ = os.environ if environ is None else environ
    user = environ.get('USER') or environ.get('LOGNAME') or 'USER'
    job = environ.get('SLURM_JOB_ID')
    if job:      # on a compute node: go through the login node
        node = environ.get('SLURMD_NODENAME') or host or socket.gethostname()
        login = environ.get('SLURM_SUBMIT_HOST') or 'LOGIN_NODE'
        where = f'{node} (Slurm job {job})'
        tunnel = f'ssh -N -J {user}@{login} -L {port}:{socket_path} {user}@{node}'
        note = (f'Use the names you normally connect with in place of {login} and {node}. If your SSH configuration '
                f'already reaches {node} through the login node, this is enough:\n'
                f'    ssh -N -L {port}:{socket_path} {node}\n'
                'The interface stops when its job ends; the runs it submitted continue.')
    else:
        host = host or socket.gethostname()
        where = host
        tunnel = f'ssh -N -L {port}:{socket_path} {user}@{host}'
        note = f'Use the name you normally connect with in place of {host}.'
    return (f'Moonstone is serving on {where}, on the private socket {socket_path}.\n'
            'No network port is open: only you can reach it, through SSH.\n\n'
            'On your own computer, run:\n'
            f'    {tunnel}\n'
            f'and open http://localhost:{port} in a browser. Leave the ssh command running; stop it with Ctrl-C.\n'
            f'{note}\n')


def build_command(argv: list) -> tuple:
    """
    (the command that serves the interface, the socket path or None, the local port of the tunnel).

    Locally the command is ``streamlit run``. With ``--remote`` it is uvicorn serving moonstone.ui.asgi on the
    Unix socket, because Streamlit's own server cannot listen on one. Other arguments go to that command.
    """
    parser = argparse.ArgumentParser(prog='moonstone-ui', add_help=False)
    parser.add_argument('--remote', action='store_true')
    parser.add_argument('--socket')
    parser.add_argument('--port', type=int, default=8501)
    options, extra = parser.parse_known_args(argv)
    if options.remote:
        socket_path = os.path.abspath(options.socket or default_socket())
        command = [sys.executable, '-m', 'uvicorn', 'moonstone.ui.asgi:app', '--uds', socket_path]
        return command + extra, socket_path, options.port
    app = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'app.py')
    command = [sys.executable, '-m', 'streamlit', 'run', app, '--server.address', '127.0.0.1',
               '--browser.gatherUsageStats', 'false',
               # one accent colour per theme, so the app still follows the system's light or dark setting
               '--theme.light.primaryColor', '#4F6FA8', '--theme.dark.primaryColor', '#8FB0E6']
    return command + extra, None, options.port


def _serve_remote(command: list, socket_path: str, port: int) -> int:
    try:
        prepare_socket(socket_path)
    except RuntimeError as error:
        print(f'Error: {error}')
        return 1
    text = connection_text(socket_path, port)
    os.makedirs(os.path.dirname(connection_file()), mode=0o700, exist_ok=True)
    with open(os.open(connection_file(), os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600), 'w') as handle:
        handle.write(text)
    print(text, flush=True)
    process = subprocess.Popen(command)
    for signum in (signal.SIGTERM, signal.SIGINT):       # stop the server with the launcher (job end, Ctrl-C)
        signal.signal(signum, lambda *_: process.terminate())
    try:
        while process.poll() is None and not os.path.exists(socket_path):
            time.sleep(0.2)
        if os.path.exists(socket_path):
            os.chmod(socket_path, 0o600)                 # uvicorn creates it open to all; the folder already protects it
        return process.wait()
    finally:
        for leftover in (socket_path, connection_file()):
            if os.path.exists(leftover):
                os.remove(leftover)


def main(argv=None) -> int:
    """Start the interface; extra arguments are passed to ``streamlit run``."""
    argv = list(sys.argv[1:] if argv is None else argv)
    if '--where' in argv:
        if os.path.exists(connection_file()):
            with open(connection_file()) as handle:
                print(handle.read())
            return 0
        print('Moonstone is not serving with --remote for this account (no connection file).')
        return 1
    try:
        import streamlit  # noqa: F401
    except ImportError:
        print('The interface needs Streamlit: pip install "moonstone[ui]"')
        return 1
    command, socket_path, port = build_command(argv)
    if socket_path is None:
        return subprocess.call(command)
    return _serve_remote(command, socket_path, port)


if __name__ == '__main__':
    sys.exit(main())
