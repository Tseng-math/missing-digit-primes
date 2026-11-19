# Missing-digit primes in bases six through eleven

**Tseng** · [Paper (PDF)](missing_digit_primes.pdf) · [LaTeX source](missing_digit_primes.tex)

For each single omission in bases 6–9, and each of the 45 pairs of omitted decimal digits, the paper obtains two-sided bounds of order `A(x)/log x`, where `A(x)` counts admissible positive integers below the real endpoint `x`. These five families share a grouped-product sieve; the available comparisons and the resulting sieve losses differ. For one omission in base 10 or 11 the result is an asymptotic equivalent, with the local factor from the allowed units digits. All counts use ordinary expansions, without leading zeros.

## Build the paper

Run `tectonic missing_digit_primes.tex`, or run `latexmk -xelatex missing_digit_primes.tex` with a TeX installation providing the Latin Modern OpenType fonts. The bibliography is embedded in the source; BibTeX, external figures and numerical verification are not required to compile it.

## Reproduce the finite calculations

From a writable extracted copy, with Python 3.12 and a C++17 compiler supporting `unsigned __int128`:

```sh
python3 -IBS supplement/verify.py all
```

The script requires only the Python standard library, with assertions enabled. It invokes `c++` by default; set `CXX` to another compiler executable if needed. No network connection is required. Exact nodal and dual witnesses are included as readable data in the Python sources. Intermediate files are created in a temporary directory and removed on completion.

To check a single group, replace `all` with `base6`, `base7`, `base8`, `base9`, `base10`, `base11`, `decimal_two`, or `common_sieve`. The full run checks the digital inequalities, rational systems and upper sums, and ends with `PASS: all numerical certificates verified.` It verifies the finite premises of the proofs, not the analytic arguments by computer.

AI Disclosure: The author selected the research topic; artificial intelligence completed all mathematical derivations, proofs, and writing.
