#pragma once

#include <cstddef>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <memory>
#include <nfl.hpp>

namespace util {
/// Helper functions for messages conversion

/**
 * Center op1 modulo op2
 * @param rop     result
 * @param op1     number op1 already reduced modulo op2, i.e. such that 0 <= op1
 * < op2-1
 * @param op2     modulus
 * @param op2Div2 floor(modulus/2)
 */
inline void center(mpz_t & rop, mpz_t const &op1, mpz_t const &op2,
		mpz_t const &op2Div2) {
	mpz_set(rop, op1);
	if (mpz_cmp(op1, op2Div2) > 0) {
		mpz_sub(rop, rop, op2);
	}
}

/* Exact comparison of ring elements. NFLlib's operator== is element-wise and
 * its operator bool is "some coefficient is non-zero", so `a == b` holds as
 * soon as a and b agree in one of the nmoduli * degree values they store,
 * which is far too weak for a verification equation. Both arguments have to be
 * in the same domain. operator!= is fine as it stands. */
template <class P> inline bool equal(P const &a, P const &b) {
	for (size_t cm = 0; cm < P::nmoduli; cm++) {
		for (size_t i = 0; i < P::degree; i++) {
			if (a(cm, i) != b(cm, i)) {
				return false;
			}
		}
	}
	return true;
}

/* Inversion in R_q. The modulus satisfies q = 1 mod 2d, so x^d + 1 splits into
 * d linear factors and a polynomial in the NTT domain already is its own CRT
 * decomposition: inverting it is inverting each slot, and it is invertible
 * exactly when no slot is zero. Both arguments are in the NTT domain, and they
 * may be the same polynomial. Returns 0 on a zero divisor, leaving inv alone.
 *
 * This replaces a dense inversion through FLINT -- convert to GMP, build a
 * polynomial mod x^d + 1, call invmod, convert back -- which cost about 95
 * million cycles and could not report a zero divisor without the caller
 * checking the return of fmpz_mod_poly_invmod, which no caller did. */
template <class P> inline int invert(P &inv, P const &a) {
	using T = typename P::value_type;
	nfl::ops::mulmod<T, nfl::simd::serial> mulmod;

	for (size_t cm = 0; cm < P::nmoduli; cm++) {
		const T q = P::get_modulus(cm);
		for (size_t i = 0; i < P::degree; i++) {
			T x = a(cm, i), r = 1;

			if (x == 0) {
				return 0;
			}
			/* Fermat: x^(q - 2) = x^-1, since q is prime. */
			for (T e = q - 2; e != 0; e >>= 1) {
				if (e & 1) {
					r = mulmod(r, x, cm);
				}
				x = mulmod(x, x, cm);
			}
			inv(cm, i) = r;
		}
	}
	return 1;
}

template <class P> inline bool is_zero(P const &a) {
	for (size_t cm = 0; cm < P::nmoduli; cm++) {
		for (size_t i = 0; i < P::degree; i++) {
			if (a(cm, i) != 0) {
				return false;
			}
		}
	}
	return true;
}
} // namespace util
