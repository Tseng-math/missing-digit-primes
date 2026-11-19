"""Build C++ certificates in a temporary directory alongside these sources."""

from contextlib import contextmanager
from os import environ
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parent


@contextmanager
def executable(source):
    with TemporaryDirectory(dir=ROOT, prefix='.verify-') as temp:
        binary = Path(temp) / 'verify'
        env = dict(environ, TMPDIR=temp, TMP=temp, TEMP=temp, CLANG_MODULE_CACHE_PATH=temp)
        run(
            [environ.get('CXX', 'c++'), '-std=c++17', '-O3', '-pipe', str(ROOT / source),
             '-o', str(binary)],
            check=True, cwd=temp, env=env,
        )
        yield binary, env
