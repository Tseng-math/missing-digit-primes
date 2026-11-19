"""Recompute the non-sieve certificates for all decimal digit pairs."""
import exponents
import moments
import projection
import pure


if __name__ == '__main__':
    moments.verify('moments', 'projection')
    projection.verify()
    exponents.verify()
    pure.verify()
