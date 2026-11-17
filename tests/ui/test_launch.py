import http.client
import os
import socket
import stat
import subprocess
import sys
import tempfile
import time
from unittest import TestCase

from moonstone.ui import launch


class UnixHTTPConnection(http.client.HTTPConnection):
    """HTTP over a Unix socket."""

    def __init__(self, path):
        super().__init__('localhost', timeout=10)
        self.socket_path = path

    def connect(self):
        self.sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        self.sock.settimeout(10)
        self.sock.connect(self.socket_path)


class TestRemoteLaunch(TestCase):

    def setUp(self):
        # a short path: Unix socket paths are limited to about 100 characters
        self.tmp = tempfile.mkdtemp(prefix='ms', dir='/tmp')
        self.folder = os.path.join(self.tmp, 'private')
        self.socket = os.path.join(self.folder, 'ui.sock')

    def test_local_command_is_unchanged(self):
        command, socket_path, port = launch.build_command(['--server.port', '8600'])
        self.assertIsNone(socket_path)
        self.assertEqual(command[command.index('--server.address') + 1], '127.0.0.1')
        self.assertEqual(command[-2:], ['--server.port', '8600'])
        self.assertNotIn('--server.headless', command)

    def test_remote_command_serves_on_the_socket_only(self):
        command, socket_path, port = launch.build_command(['--remote', '--socket', self.socket, '--port', '8600'])
        self.assertEqual((socket_path, port), (self.socket, 8600))
        self.assertEqual(command[1:5], ['-m', 'uvicorn', 'moonstone.ui.asgi:app', '--uds'])
        self.assertEqual(command[5], self.socket)
        for network_option in ('127.0.0.1', '--host', '--port', '--server.address'):
            self.assertNotIn(network_option, command)

    def test_default_socket_is_in_a_private_folder_of_the_home(self):
        self.assertEqual(os.path.dirname(launch.default_socket()), os.path.expanduser('~/.moonstone'))
        self.assertIn(socket.gethostname(), os.path.basename(launch.default_socket()))

    def test_folder_is_created_private(self):
        launch.prepare_socket(self.socket)
        self.assertEqual(stat.S_IMODE(os.stat(self.folder).st_mode), 0o700)

    def test_folder_with_wider_permissions_is_refused(self):
        os.makedirs(self.folder, mode=0o755)
        os.chmod(self.folder, 0o755)
        with self.assertRaisesRegex(RuntimeError, 'other users'):
            launch.prepare_socket(self.socket)

    def test_stale_socket_is_replaced_and_live_socket_is_an_error(self):
        launch.prepare_socket(self.socket)
        server = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        server.bind(self.socket)
        server.listen(1)
        with self.assertRaisesRegex(RuntimeError, 'already serving'):
            launch.prepare_socket(self.socket)
        server.close()                               # the file stays: a stale socket
        self.assertTrue(os.path.exists(self.socket))
        launch.prepare_socket(self.socket)
        self.assertFalse(os.path.exists(self.socket))

    def test_connection_instructions(self):
        text = launch.connection_text('/home/u/.moonstone/ui.sock', 8501, {'USER': 'ada'}, host='cluster-login')
        self.assertIn('ssh -N -L 8501:/home/u/.moonstone/ui.sock ada@cluster-login', text)
        self.assertIn('http://localhost:8501', text)
        self.assertNotIn(' -J ', text)
        in_job = launch.connection_text('/home/u/.moonstone/ui.sock', 8600, {
            'USER': 'ada', 'SLURM_JOB_ID': '42', 'SLURMD_NODENAME': 'node-7', 'SLURM_SUBMIT_HOST': 'cluster-login'})
        self.assertIn('ssh -N -J ada@cluster-login -L 8600:/home/u/.moonstone/ui.sock ada@node-7', in_job)
        self.assertIn('Slurm job 42', in_job)

    def test_where_without_a_server(self):
        result = subprocess.run([sys.executable, '-m', 'moonstone.ui.launch', '--where'], capture_output=True, text=True,
                                env=dict(os.environ, HOME=self.tmp))
        self.assertEqual(result.returncode, 1)
        self.assertIn('not serving', result.stdout)

    def test_remote_serves_on_a_private_socket(self):
        """The real server: reachable through the socket, described in connection.txt, cleaned up when it stops."""
        process = subprocess.Popen([sys.executable, '-m', 'moonstone.ui.launch', '--remote', '--socket', self.socket],
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                                   env=dict(os.environ, HOME=self.tmp))
        try:
            status = None
            for _ in range(90):
                try:
                    connection = UnixHTTPConnection(self.socket)
                    connection.request('GET', '/_stcore/health')
                    status = connection.getresponse().status
                    break
                except OSError:
                    self.assertIsNone(process.poll(), 'the server stopped')
                    time.sleep(1)
            self.assertEqual(status, 200)
            self.assertEqual(stat.S_IMODE(os.stat(self.folder).st_mode), 0o700)      # what keeps other users out
            for _ in range(20):
                if stat.S_IMODE(os.stat(self.socket).st_mode) == 0o600:
                    break
                time.sleep(0.5)
            self.assertEqual(stat.S_IMODE(os.stat(self.socket).st_mode), 0o600)
            connection = UnixHTTPConnection(self.socket)
            connection.request('GET', '/')
            self.assertIn(b'<title>', connection.getresponse().read())               # the interface itself
            where = subprocess.run([sys.executable, '-m', 'moonstone.ui.launch', '--where'], capture_output=True,
                                   text=True, env=dict(os.environ, HOME=self.tmp))
            self.assertIn(f'-L 8501:{self.socket}', where.stdout)
            try:
                import psutil
                listening = [c for child in [psutil.Process(process.pid)] + psutil.Process(process.pid).children(True)
                             for c in child.net_connections(kind='tcp') if c.status == psutil.CONN_LISTEN]
                self.assertEqual(listening, [])                                      # no TCP port at all
            except ImportError:
                pass
        finally:
            process.terminate()
            output = process.communicate(timeout=30)[0]
        self.assertIn('ssh -N', output)
        self.assertFalse(os.path.exists(self.socket))
        self.assertFalse(os.path.exists(os.path.join(self.tmp, '.moonstone', 'connection.txt')))
