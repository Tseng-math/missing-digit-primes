"""Run all finite checks, or select a group. Use python3 -IBS supplement/verify.py --help."""
from argparse import ArgumentParser
import ast
from hashlib import sha256
import os
from pathlib import Path
from shutil import copyfile
from subprocess import run
import sys
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parent
CHECKS = {
    'base6': ('cert/verify_moments.py', 'cert/verify_mixed.py', 'cert/verify_pure_refined.py',
              'bilinear/verify_projection.py', 'bilinear/verify_pure_tail.py', 'adaptive/verify_gram.py',
              'small_periods/verify_small_period.py', 'positive/verify_mask.py', 'check_weights.py'),
    'base7': ('verify_moments.py', 'verify_orbits.py', 'verify_endpoints.py'),
    'base8': ('verify_first_moment.py', 'verify_near_parseval.py', 'verify_exponents.py', 'verify_loss.py'),
    'base9': ('verify_analytic.py',),
    'base10': ('verify_base10.py',),
    'base11': ('verify_base11.py',),
    'decimal_two': ('verify_analytic.py',),
}
WINDOW_FOR = {'base6': 'six', 'base7': 'seven_pairs', 'decimal_two': 'seven_pairs', 'base9': 'nine'}
WITNESS_HASHES = {'base10/forward_certificates.json': '6d3f515d38d2879aeb50b1132dea3dd188693f017858b537e8ed0c6d79b7053c',
 'base10/hellinger_certificates.json': 'cd77f39924e3d81d98402d6b97e88282b48b20c8a833b5de45c19592b365f367',
 'base10/moment_classical_d0.json': '481f8314c3ec1227637fb0d67facb3550a407409eaf2e72174e52ee3c56d3ef5',
 'base10/moment_classical_d1.json': 'eb6e951888c400c61a6b00b5d1445a34af24324afc78aca472870934368170c8',
 'base10/moment_classical_d2.json': '19a522a396a7e893aeb6a36b66ba2c997e0120233bdffa41089dd4966f1b79bc',
 'base10/moment_classical_d3.json': '8d6c4e5278a25613a5c9c05936dfdd6c366d2894b2155756a5bb982b06f0f15e',
 'base10/moment_classical_d4.json': 'd1246142fd3cddbc07c3903567abdbe5a6a72ada96cc3170b658285f63453a2c',
 'base10/moment_one_d0.json': 'c033f3074be15bf817b09d9f39c1908e459853af4620886390f9feb517dbdacf',
 'base10/moment_one_d1.json': 'ee0f68504064ad7f8ed86adb9ddf7a0c477a3766e27114b302d90ba51540c139',
 'base10/moment_one_d2.json': '3ecdc0a984a6c53a17f33b6b12243bb997719b0980f692f3bb2138867a6378e8',
 'base10/moment_one_d3.json': 'e90ef289325ca8612edce860a94fe27fa5133ce38fdd48866001f8e8384ff08d',
 'base10/moment_one_d4.json': 'cd0fd1e106b1fc2cc96b72040228d3a266a7796e995bdfc571c917e8ff814f47',
 'base10/moment_three_halves_d0.json': '9b896de03de7d1eef75726f7dd4c21e6a9eaf226ff398c6f226884c1a48edfcf',
 'base10/moment_three_halves_d1.json': 'ce8e17f5aefe54311f28696076e4d260397871d1151000d7894b7d7b60accc97',
 'base10/moment_three_halves_d2.json': '22fae81cd76eec03924fa73675fe2c4ccaee2765dcca63c6b71b81703e137fc1',
 'base10/moment_three_halves_d3.json': '9ae8aae413e6e45eecafa35d0d276694a9aba85775674bd91cad958c1b2a8c84',
 'base10/moment_three_halves_d4.json': '06a6b316a3884bf1473f45f7611b1450fbed683fce9d4601353acd5d198e3705',
 'base10/reverse_central_certificates.json': 'a183cb6910e82f79e5ec634c689778f44c15d592699fc487568892d13912f9b8',
 'base10/reverse_uniform_certificates.json': '88749d6d8f7225fd84537a70cf400ae000791b08abad5ec61bdfb6bd6a6f83e7',
 'base11/moments.json': '31999f0986fdf0b84b0ea996aba0b16b82e8af7069d807f2c29fd52e13880df8',
 'base8/first_moment_nodes.json': 'd7eaa7bcc26111c287999192cfaf06170d5fe09d28878d28c4cda23b663ae921',
 'base8/near_parseval_nodes.json': '63ff36e59ea437f5789c5ea8b2cde02c911a440286b5bdbc1eaca15421aef52a',
 'base9/moment_certificates.json': '115c2fe032de0e780d2ce31de71e2cbc5cde7143e90bf7adf8640adfcaae7822'}


def prepare(work):
    restored = {}
    for group in (*CHECKS, 'common_sieve'):
        for source in sorted((ROOT/group).rglob('*')):
            if not source.is_file():
                continue
            if source.suffix not in ('.py', '.cpp'):
                raise SystemExit(f'Unexpected file in supplement: {source}')
            dest = work/source.relative_to(ROOT)
            dest.parent.mkdir(parents=True, exist_ok=True)
            if source.name != 'witnesses.py':
                copyfile(source, dest)
                continue
            assignment, = ast.parse(source.read_text(encoding='utf-8')).body
            for name, text in ast.literal_eval(assignment.value).items():
                path = Path(group)/name
                data = text.encode('utf-8')
                if WITNESS_HASHES.get(str(path)) != sha256(data).hexdigest():
                    raise SystemExit(f'Changed or unexpected witness: {path}')
                (work/path).write_bytes(data)
                restored[str(path)] = True
    if restored.keys() != WITNESS_HASHES.keys():
        raise SystemExit('Incomplete witness set')
    print(f'Restored and checked {len(restored)} transient witness files.', flush=True)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('group', nargs='?', default='all', choices=('all', *CHECKS, 'common_sieve'))
    parser.add_argument('--inputs-only', action='store_true')
    args = parser.parse_args()
    if not __debug__:
        parser.error('Run with assertions enabled.')
    with TemporaryDirectory(dir=ROOT.parent, prefix='.certificate-') as directory:
        work = Path(directory)
        prepare(work)
        if args.inputs_only:
            return
        env = {key:value for key,value in os.environ.items() if not key.startswith('PYTHON')}
        env.update(TMPDIR=directory, TMP=directory, TEMP=directory, CLANG_MODULE_CACHE_PATH=directory)

        def check(program, *arguments):
            print(f'Running {program} {" ".join(arguments)}', flush=True)
            run([sys.executable, '-EBuS', str(work/program), *arguments], cwd=work, env=env, check=True)

        groups = tuple(CHECKS) if args.group == 'all' else (args.group,)
        windows = tuple(dict.fromkeys(WINDOW_FOR.values())) if args.group == 'common_sieve' else tuple(
            dict.fromkeys(WINDOW_FOR[group] for group in groups if group in WINDOW_FOR))
        if windows:
            check('common_sieve/verify.py', *windows)
        for group in groups:
            for script in CHECKS.get(group, ()):
                check(f'{group}/{script}')
            print(f'PASS {group}', flush=True)
    print(f'PASS: {args.group} numerical certificates verified.', flush=True)


if __name__ == '__main__':
    main()
