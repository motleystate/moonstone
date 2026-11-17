.. _remote:

**********************************
Using the interface on a cluster
**********************************

Reads and databases live on the cluster, so the interface runs there too, and you open it in the browser of your own
computer through an SSH tunnel. Nothing is installed on your computer: it needs only ``ssh`` and a browser.

Why a private socket and not a port
===================================

On a machine shared with other users, a web server listening on a network port can be reached by any of them, even
when it is bound to ``127.0.0.1``. The interface can start jobs and read files as you, so it must not be reachable
by anyone else.

``moonstone-ui --remote`` therefore opens **no network port**. It listens on a Unix socket, a file, placed in a
folder that only you can enter (``~/.moonstone``, mode 700). SSH connects a port of *your* computer to that file.

Starting it
===========

On the cluster, in the environment where Moonstone is installed:

.. code-block:: bash

    moonstone-ui --remote

It prints what to run on your computer, and writes the same text to ``~/.moonstone/connection.txt``:

.. code-block:: text

    Moonstone is serving on node-12 (Slurm job 123456), on the private socket /home/ada/.moonstone/ui-node-12.sock.
    No network port is open: only you can reach it, through SSH.

    On your own computer, run:
        ssh -N -J ada@login-node -L 8501:/home/ada/.moonstone/ui-node-12.sock ada@node-12
    and open http://localhost:8501 in a browser. Leave the ssh command running; stop it with Ctrl-C.

``moonstone-ui --where`` prints that text again, from any node.

Where to start it depends on the rules of your cluster:

In a Slurm job (preferred)
    The interface runs on a compute node, like any other work:

    .. code-block:: bash

        sbatch -c 2 --mem 8G --wrap "moonstone-ui --remote"
        moonstone-ui --where          # once the job has started

    This needs SSH access from the login node to a node where you have a job, which most clusters allow.

On the login node
    Only if your cluster allows a long-lived process there. The command for your computer is then simply
    ``ssh -N -L 8501:SOCKET USER@LOGIN_NODE``.

Start it from the folder where the run folders should be created (``moonstone_runs`` by default), or set the run
directory in the sidebar.

On your computer
================

Run the ``ssh`` command that was printed, and open ``http://localhost:8501``. The command is the same on Linux, macOS
and Windows (OpenSSH is included in Windows 10 and later). If port 8501 is taken on your computer, start the
interface with ``--port`` (for example ``moonstone-ui --remote --port 8600``) or change the first number of ``-L``.

With an SSH configuration that already reaches the compute node through the login node, the command is shorter:

.. code-block:: text

    # ~/.ssh/config on your computer
    Host mycluster
        HostName login.mycluster.example
        User ada

    Host node-*
        User ada
        ProxyJump mycluster

.. code-block:: bash

    ssh -N -L 8501:/home/ada/.moonstone/ui-node-12.sock node-12

When things stop
================

* **You close the browser or the tunnel:** nothing is lost. Runs submitted to Slurm continue; reconnect and find them
  under *Previous runs*, with their current state.
* **The job of the interface ends:** the same. Start the interface again; the runs are read from their folders.
* **A run submitted without "Submit to Slurm"** runs inside the interface and stops with it. On a cluster, tick it.

If something does not work
==========================

``Moonstone is already serving on ...``
    An interface of yours is running on that machine: ``moonstone-ui --where`` says how to reach it.

``... can be entered by other users``
    The folder of the socket is not private. Run the ``chmod 700`` command that the message gives.

``REMOTE HOST IDENTIFICATION HAS CHANGED`` from ``ssh``
    The node presents another key than the one your computer remembers. Nodes are sometimes reinstalled, but do not
    ignore the message: check with the administrators of the cluster before accepting the new key.

The page loads but stays grey
    The tunnel was closed. Run the ``ssh`` command again and reload the page.
