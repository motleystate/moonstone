*************************
Contributing to Moonstone
*************************

This page is under construction and will soon guide you to contribute to the project.

.. Note::
    In the meantime, for any request or suggestion, you can refer to the issues_ on the github repository.

.. _issues: https://github.com/motleystate/moonstone/issues

Releasing a new version
=======================

The online documentation (Read the Docs) is built from the version of Moonstone published on the Python Package Index (PyPI), because ``docs/requirements.txt`` installs ``moonstone`` from PyPI. It therefore describes the latest release, not unreleased changes on ``master``. Release and documentation must be updated together:

1. Update the version in ``moonstone/version.py`` and ``pyproject.toml``, and add the changelog entry in ``docs/changelogs.rst``.
2. Create a GitHub release with a matching tag. The ``Upload Python Package`` workflow publishes it to PyPI.
3. Check that the new version is on PyPI before going further.
4. If the minimum Python version changed, update ``python`` (and ``os`` if needed) in ``.readthedocs.yaml`` to a version the new release supports.
5. On the Read the Docs project, make sure the default version is ``stable`` (the newest tag), then check the build log and the online pages for the new version.
