"""Run the exact verification programs for bases eight through eleven."""

from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
CHECKS = (
    'base8/verify_first_moment.py',
    'base8/verify_near_parseval.py',
    'base8/verify_exponents.py',
    'base9/verify_base9.py',
    'base10/verify_base10.py',
    'base11/verify_base11.py',
)


def main():
    if sys.flags.optimize:
        raise SystemExit('Run the verification programs with assertions enabled.')
    for name in CHECKS:
        print(f'Verifying {name}', flush=True)
        subprocess.run([sys.executable, '-B', '-S', str(ROOT/name)],
                       cwd=ROOT, check=True)
    print('All numerical certificates verified.')


if __name__ == '__main__':
    main()
