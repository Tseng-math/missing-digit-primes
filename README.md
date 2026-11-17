# Missing-digit primes in bases eight through eleven

**Tseng**

The paper proves two-sided prime counts of the expected order of magnitude in bases 8 and 9, and exact asymptotic formulas in bases 10 and 11. Both results hold for every omitted digit and arbitrary real endpoints.

- [Manuscript (PDF)](missing_digit_primes_8_11.pdf)
- [Complete LaTeX source](missing_digit_primes_8_11.tex)
- [Verification programs and certificate data](supplement/)

## Build the PDF

With Tectonic installed, run the following from the repository root. Tectonic may download its TeX bundle on the first build.

```sh
mkdir -p build
tectonic --untrusted --outdir build missing_digit_primes_8_11.tex
```

The rebuilt PDF is written to `build/missing_digit_primes_8_11.pdf`.

## Reproduce the certificate checks

The verifiers use Python 3.8 or later and its standard library. The base-nine check also requires `g++` with C++17 and `unsigned __int128` support. Run with assertions enabled from a writable copy of the repository.

```sh
python3 -B -S supplement/verify.py
```

The command recomputes the finite certificates for bases 8 through 11 and prints `All numerical certificates verified.` on success. The base-nine verifier removes its temporary build files after the check. The individual scripts can also be run with `python3 -B -S` followed by the script path.

| Base | Certificate checks |
| --- | --- |
| 8 | [First moment](supplement/base8/verify_first_moment.py), [moment of order 127/64](supplement/base8/verify_near_parseval.py), and [exponent and sieve inequalities](supplement/base8/verify_exponents.py) |
| 9 | [Fourier moments and complete-tuple sieve loss](supplement/base9/verify_base9.py) |
| 10 | [Moments, affine correlations, and Type II exponent inequalities](supplement/base10/verify_base10.py) |
| 11 | [Moments and affine correlations](supplement/base11/verify_base11.py) |

## Citation

Tseng. *Missing-digit primes in bases eight through eleven*. Manuscript.

```bibtex
@unpublished{TsengMissingDigitPrimes,
  author = {Tseng},
  title  = {Missing-digit primes in bases eight through eleven},
  note   = {Manuscript}
}
```

Machine-readable citation metadata is available in [CITATION.cff](CITATION.cff).

## AI Disclosure

The author selected the research topic; artificial intelligence completed all mathematical derivations, proofs, and writing.
