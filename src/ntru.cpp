#include "test.h"
#include "bench.h"
#include "common.h"
#include <sys/random.h>
#include <math.h>
#include <stdlib.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpz_mod_poly.h>
#include "blake3.h"
#include "assert.h"
#include "sample_z_small.h"


void ntru_sample_message(ntru_params::poly_p &m) {
    // Same as BGV, but we are going to use different params so created a specific one.
    std::array<mpz_t, NTRU_DEGREE> coeffs;
    uint64_t buf;
    size_t bits_in_moduli_product = ntru_params::poly_p::bits_in_moduli_product();
    for (size_t i = 0; i < ntru_params::poly_p::degree; i++) {
        mpz_init2(coeffs[i], bits_in_moduli_product << 2);
    }
    for (size_t j = 0; j < ntru_params::poly_p::degree; j += 32) {
        ssize_t grret = getrandom(&buf, sizeof(buf), 0);
        if (grret == -1) {
            throw "Could not generate random seed. Check if something is wrong with the system prng.";
        }
        for (size_t k = 0; k < 64; k += 2) {
            mpz_set_ui(coeffs[j + k / 2], (buf >> k) % NTRU_PRIMEP);
        }
    }
    m.mpz2poly(coeffs);
    for (size_t i = 0; i < ntru_params::poly_p::degree; i++) {
        mpz_clear(coeffs[i]);
    }
}

/**
 * @brief Computes the multiplicative inverse of a polynomial in a given field.
 *
 * @param inv Output parameter to store the computed inverse polynomial.
 * @param p Input polynomial for which the inverse is to be computed.
 */
void poly_inverse(ntru_params::poly_q &inv, ntru_params::poly_q p) {
    // Declare an array to store coefficients of the polynomial
    std::array<mpz_t, ntru_params::poly_q::degree> coeffs;
    // Declare a variable to store the modulus of the field
    fmpz_t q;
    // Declare variables to store the polynomial and the irreducible polynomial
    fmpz_mod_poly_t poly, irred;
    // Declare a context variable for modular arithmetic operations
    fmpz_mod_ctx_t ctx_q;

    // Initialize the modulus variable
    fmpz_init(q);
    // Initialize the coefficients array with the appropriate bit size
    for (size_t i = 0; i < ntru_params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], (ntru_params::poly_q::bits_in_moduli_product() << 2));
    }

    // Set the modulus value to the product of moduli for the polynomial
    fmpz_set_mpz(q, ntru_params::poly_q::moduli_product());
    // Initialize the context for modular arithmetic with the modulus
    fmpz_mod_ctx_init(ctx_q, q);
    // Initialize the polynomial and irreducible polynomial variables in the context
    fmpz_mod_poly_init(poly, ctx_q);
    fmpz_mod_poly_init(irred, ctx_q);

    // Convert the polynomial to its coefficient representation
    p.poly2mpz(coeffs);
    // Define the irreducible polynomial for the field
    fmpz_mod_poly_set_coeff_ui(irred, ntru_params::poly_q::degree, 1, ctx_q);
    fmpz_mod_poly_set_coeff_ui(irred, 0, 1, ctx_q);

    // Set the polynomial coefficients from the array
    for (size_t i = 0; i < ntru_params::poly_q::degree; i++) {
        fmpz_mod_poly_set_coeff_mpz(poly, i, coeffs[i], ctx_q);
    }
    // Compute the multiplicative inverse of the polynomial modulo the irreducible polynomial
    fmpz_mod_poly_invmod(poly, poly, irred, ctx_q);

    // Retrieve the coefficients of the inverse polynomial
    for (size_t i = 0; i < ntru_params::poly_q::degree; i++) {
        fmpz_mod_poly_get_coeff_mpz(coeffs[i], poly, i, ctx_q);
    }

    // Convert the coefficient representation back to the polynomial form
    inv.mpz2poly(coeffs);

    // Clear the memory allocated for the modulus
    fmpz_clear(q);
    // Clear the memory allocated for the coefficients array
    for (size_t i = 0; i < ntru_params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
}

/*
 * Test norm of polynomial less than some bound sqrRoot(t^2*sigma^2*degree)
 * */
bool ntru_test_norm(ntru_params::poly_q r, double sigma_sqr, double t) {
    array<mpz_t, ntru_params::poly_q::degree> coeffs;
    mpz_t norm, qDivBy2, tmp;
    mpz_inits(norm, qDivBy2, tmp, nullptr);
    for (size_t i = 0; i < ntru_params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], ntru_params::poly_q::bits_in_moduli_product() << 2);
    }
    r.poly2mpz(coeffs);
    mpz_fdiv_q_2exp(qDivBy2, ntru_params::poly_q::moduli_product(), 1);
    mpz_set_ui(norm, 0);
    for (size_t i = 0; i < ntru_params::poly_q::degree; i++) {
        util::center(coeffs[i], coeffs[i], ntru_params::poly_q::moduli_product(), qDivBy2);
        mpz_mul(tmp, coeffs[i], coeffs[i]);
        mpz_add(norm, norm, tmp);
    }
    long double bound = t * t * sigma_sqr * ntru_params::poly_q::degree;
    int result = mpz_cmp_ui(norm, bound) < 0;
    mpz_clears(norm, qDivBy2, tmp, nullptr);
    for (size_t i = 0; i < ntru_params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
    return result;
}

void ntru_keygen(ntru_params::poly_q &pk, ntru_params::poly_q &sk) {
    ntru_params::poly_q f, g, f_inv;
    array<mpz_t, ntru_params::poly_q::degree> coeffs_f, coeffs_g, coeffs;

    for (size_t i = 0; i < ntru_params::poly_q::degree; i++) {
        mpz_init2(coeffs_f[i], ntru_params::poly_q::bits_in_moduli_product() << 2);
        mpz_init2(coeffs_g[i], ntru_params::poly_q::bits_in_moduli_product() << 2);
        mpz_init2(coeffs[i], ntru_params::poly_q::bits_in_moduli_product() << 2);
    }

    do {
        for (size_t k = 0; k < ntru_params::poly_q::degree; k++) {
            int64_t coeff_f;
            do {
                coeff_f = sample_z(0.0, NTRU_SIGMA);
            } while ((k == 0 && coeff_f % NTRU_PRIMEP != 1) || (k >= 1 && coeff_f % NTRU_PRIMEP != 0));
            int64_t coeff_g = sample_z(0.0, NTRU_SIGMA);
            mpz_set_si(coeffs_f[k], coeff_f);
            mpz_set_si(coeffs_g[k], coeff_g);
        }
        f.mpz2poly(coeffs_f);
        g.mpz2poly(coeffs_g);
    } while (!ntru_test_norm(f, NTRU_SIGMA*NTRU_SIGMA, 1.058));

    poly_inverse(f_inv, f);

    f.ntt_pow_phi();
    g.ntt_pow_phi();
    f_inv.ntt_pow_phi();

    pk = g * f_inv;
    sk = f;

    for (size_t i = 0; i < ntru_params::poly_q::degree; i++) {
        mpz_clear(coeffs_f[i]);
        mpz_clear(coeffs_g[i]);
        mpz_clear(coeffs[i]);
    }
}

void ntru_keyshare(ntru_params::poly_q s[], size_t shares, ntru_params::poly_q &sk) {
    ntru_params::poly_q t = sk;
    for (size_t i = 1; i < shares; i++) {
        s[i] = nfl::uniform();
        s[i].ntt_pow_phi();
        t = t - s[i];
    }
    s[0] = t;
}

void ntru_encrypt(ntru_params::poly_q &c, ntru_params::poly_q &pk, ntru_params::poly_p &m) {
    ntru_params::poly_q s = nfl::ZO_dist();
    ntru_params::poly_q e = nfl::ZO_dist();
    ntru_params::poly_q m_;
    std::array<mpz_t, NTRU_DEGREE> coeffs;

    for (int i = 0; i < ntru_params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], ntru_params::poly_q::bits_in_moduli_product() << 2);
    }

    s.ntt_pow_phi();
    e.ntt_pow_phi();

    m.poly2mpz(coeffs);
    m_.mpz2poly(coeffs);

    m_.ntt_pow_phi();

    c = (pk * s + e) + (pk * s + e) + m_;

    for (int i = 0; i < ntru_params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
}

void ntru_add(ntru_params::poly_q &c, ntru_params::poly_q &c1, ntru_params::poly_q &c2) {
    c = c1 + c2;
}

void ntru_decrypt(ntru_params::poly_p &m, ntru_params::poly_q &c, ntru_params::poly_q &sk) {
    std::array<mpz_t, NTRU_DEGREE> coeffs;
    ntru_params::poly_q t = sk * c;
    mpz_t qDivBy2;

    mpz_init(qDivBy2);
    for (size_t i = 0; i < ntru_params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], ntru_params::poly_q::bits_in_moduli_product() << 2);
    }

    mpz_fdiv_q_2exp(qDivBy2, ntru_params::poly_q::moduli_product(), 1);
    t.invntt_pow_invphi();
    t.poly2mpz(coeffs);
    for (size_t i = 0; i < ntru_params::poly_q::degree; i++) {
        util::center(coeffs[i], coeffs[i], ntru_params::poly_q::moduli_product(), qDivBy2);
        mpz_mod_ui(coeffs[i], coeffs[i], NTRU_PRIMEP);
    }
    m.mpz2poly(coeffs);

    mpz_clear(qDivBy2);
    for (size_t i = 0; i < ntru_params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
}

void ntru_distdec(ntru_params::poly_q &dsj, ntru_params::poly_q &ci, ntru_params::poly_q &dkj) {
    std::array<mpz_t, NTRU_DEGREE> coeffs;
    ntru_params::poly_q dsij, Ej;
    mpz_t qDivBy2, bound;

    mpz_init(qDivBy2);
    mpz_init(bound);
    for (size_t i = 0; i < ntru_params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], ntru_params::poly_q::bits_in_moduli_product() << 2);
    }

    mpz_fdiv_q_2exp(qDivBy2, ntru_params::poly_q::moduli_product(), 1);
    mpz_set_str(bound, NTRU_BOUND_D, 10);

    Ej = nfl::uniform();
    Ej.poly2mpz(coeffs);
    for (size_t i = 0; i < ntru_params::poly_q::degree; i++) {
        util::center(coeffs[i], coeffs[i], ntru_params::poly_q::moduli_product(), qDivBy2);
        mpz_mod(coeffs[i], coeffs[i], bound);
    }
    Ej.mpz2poly(coeffs);
    Ej.ntt_pow_phi();
    dsij = dkj * ci;
    dsj = dsij;
    for (size_t i = 0; i < NTRU_PRIMEP; i++) {
        dsj = dsj + Ej;
    }
}

void ntru_comb(ntru_params::poly_p &m, ntru_params::poly_q c, ntru_params::poly_q dsj[], size_t shares) {
    std::array<mpz_t, NTRU_DEGREE> coeffs;
    ntru_params::poly_q v;
    mpz_t qDivBy2;

    mpz_init(qDivBy2);
    for (size_t i = 0; i < ntru_params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], ntru_params::poly_q::bits_in_moduli_product() << 2);
    }

    mpz_fdiv_q_2exp(qDivBy2, ntru_params::poly_q::moduli_product(), 1);
    v = dsj[0];
    for (size_t i = 1; i < shares; i++) {
        v = v + dsj[i];
    }
    v.invntt_pow_invphi();
    v.poly2mpz(coeffs);
    for (size_t i = 0; i < ntru_params::poly_q::degree; i++) {
        util::center(coeffs[i], coeffs[i], ntru_params::poly_q::moduli_product(), qDivBy2);
        mpz_mod_ui(coeffs[i], coeffs[i], NTRU_PRIMEP);
    }
    m.mpz2poly(coeffs);

    mpz_clear(qDivBy2);
    for (size_t i = 0; i < ntru_params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
}

#ifdef MAIN

static void ntru_test() {

    ntru_params::poly_q sk, pk, c1, c2, s[PARTIES], t[PARTIES], acc;
    ntru_params::poly_p m, _m;

    ntru_keygen(pk, sk);
    ntru_sample_message(m);

    // Test that NTRU encryption is consistent.
    TEST_BEGIN("NTRU encryption is consistent") {
        ntru_encrypt(c1, pk, m);  // Encrypt the message.
        ntru_decrypt(_m, c1, sk);  // Decrypt the ciphertext.
        TEST_ASSERT(m - _m == 0, end);  // Check that the decrypted message matches the original.

        ntru_sample_message(m);  // Sample another message.
        ntru_decrypt(_m, c1, sk);  // Decrypt the previous ciphertext.
        TEST_ASSERT(m - _m != 0, end);  // Check that the decrypted message does not match the new message.

        ntru_encrypt(c1, pk, m);  // Encrypt the new message.
        ntru_keygen(pk, sk);  // Generate a new key pair.
        ntru_decrypt(_m, c1, sk);  // Decrypt the ciphertext with the new secret key.
        TEST_ASSERT(m - _m != 0, end);  // Check that the decryption is not successful with the wrong key.
    }
    TEST_END;
    TEST_BEGIN("NTRU distributed decryption is consistent") {
        ntru_keygen(pk, sk);
        ntru_encrypt(c1, pk, m);
        ntru_keyshare(s, NTRU_PARTIES, sk);
        acc = s[0];
        for (size_t j = 1; j < NTRU_PARTIES; j++) {
            acc = acc + s[j];
        }
        TEST_ASSERT(sk - acc == 0, end);
        for (size_t j = 0; j < NTRU_PARTIES; j++) {
            ntru_distdec(t[j], c1, s[j]);
        }
        ntru_comb(_m, c1, t, NTRU_PARTIES);
        TEST_ASSERT(m - _m == 0, end);
    }
    TEST_END;

    end:
    return;
}

static void bench() {
    ntru_params::poly_q pk, sk, s[NTRU_PARTIES], t[NTRU_PARTIES], acc;
    ntru_params::poly_p m, _m;
    ntru_params::poly_q c;

    BENCH_SMALL("ntru_keygen", ntru_keygen(pk, sk));

    BENCH_BEGIN("ntru_sample_message")
        {
            BENCH_ADD(ntru_sample_message(m));
        }
    BENCH_END;

    BENCH_BEGIN("ntru_encrypt")
        {
            BENCH_ADD(ntru_encrypt(c, pk, m));
        }
    BENCH_END;

    BENCH_BEGIN("ntru_decrypt")
        {
            ntru_encrypt(c, pk, m);
            BENCH_ADD(ntru_decrypt(_m, c, sk));
        }
    BENCH_END;

    ntru_keyshare(t, NTRU_PARTIES, sk);
    BENCH_BEGIN("ntru_distdec")
        {
            ntru_encrypt(c, pk, m);
            BENCH_ADD(ntru_distdec(t[0], c, s[0]));
        }
    BENCH_END;
    for (size_t i = 1; i < PARTIES; i++) {
        ntru_distdec(t[i], c, s[i]);
    }

    BENCH_BEGIN("ntru_comb")
        {
            BENCH_ADD(ntru_comb(_m, c, t, PARTIES));
        }
    BENCH_END;
}

// Main function to run tests and benchmarks.
int main(int argc, char *arv[]) {
    printf("\n** Tests for NTRU encryption:\n\n");
    ntru_test();

    printf("\n** Benchmarks for NTRU encryption:\n\n");
    bench();

    return 0;
}

#endif
