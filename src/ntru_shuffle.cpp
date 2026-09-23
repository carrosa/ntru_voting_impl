#include <math.h>
#include <stdlib.h>

#include "flint_util.h"

#include "blake3.h"
#include "common.h"
#include "test.h"
#include "bench.h"
#include "assert.h"
#include "sample_z_small.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

/* Number of messages in the shuffle. Overridable, like TAU in ntru_pismall:
 * the arrays below are MSGS long and the proof holds all of them in memory. */
#ifndef MSGS
#define MSGS         1000
#endif

/* An params::poly_q is 16 KiB, so anything dimensioned by MSGS is far too
 * large to be a local variable: at MSGS = 1000 the buffers below add up to
 * hundreds of megabytes, which gcc reported as stack frames of 295 MB for run()
 * and 328 MB for bench(). They are allocated on the heap by shuffle_alloc() and
 * the pointers index exactly like the arrays they replace; the functions that
 * take them as arguments keep their signatures, since an array parameter is a
 * pointer either way. inv_tmp is simul_inverse's scratch space and _r, theta
 * and inv are the prover's, all equally oversized. */
static commit_t *com, *d, *cs;
static vector<params::poly_q> *r, *_r;
static params::poly_q *ms, *_ms, *s, *theta, *inv, *inv_tmp;
static params::poly_q (*y)[NTRU_WIDTH], (*_y)[NTRU_WIDTH];
static params::poly_q *t, *_t, *u;

static void shuffle_alloc(void) {
    com = new commit_t[MSGS];
    d = new commit_t[MSGS];
    cs = new commit_t[MSGS];
    r = new vector<params::poly_q>[MSGS];
    _r = new vector<params::poly_q>[MSGS];
    ms = new params::poly_q[MSGS];
    _ms = new params::poly_q[MSGS];
    s = new params::poly_q[MSGS];
    theta = new params::poly_q[MSGS];
    inv = new params::poly_q[MSGS];
    inv_tmp = new params::poly_q[MSGS];
    y = new params::poly_q[MSGS][NTRU_WIDTH];
    _y = new params::poly_q[MSGS][NTRU_WIDTH];
    t = new params::poly_q[MSGS];
    _t = new params::poly_q[MSGS];
    u = new params::poly_q[MSGS];
}

static void shuffle_free(void) {
    delete[] com;
    delete[] d;
    delete[] cs;
    delete[] r;
    delete[] _r;
    delete[] ms;
    delete[] _ms;
    delete[] s;
    delete[] theta;
    delete[] inv;
    delete[] inv_tmp;
    delete[] y;
    delete[] _y;
    delete[] t;
    delete[] _t;
    delete[] u;
}


static void lin_hash(params::poly_q &beta, comkey_t &key, commit_t x,
                     commit_t y, params::poly_q alpha[2], params::poly_q &u,
                     params::poly_q t, params::poly_q _t) {
    uint8_t hash[BLAKE3_OUT_LEN];
    blake3_hasher hasher;

    blake3_hasher_init(&hasher);

    /* Hash public key. */
    for (size_t i = 0; i < NTRU_HEIGHT; i++) {
        for (int j = 0; j < NTRU_WIDTH - NTRU_HEIGHT; j++) {
            blake3_hasher_update(&hasher, (const uint8_t *) key.A1[i][j].data(),
                                 8 * NTRU_DEGREE);
        }
    }
    for (size_t j = 0; j < NTRU_WIDTH; j++) {
        blake3_hasher_update(&hasher, (const uint8_t *) key.A2[0][j].data(),
                             8 * NTRU_DEGREE);
    }

    /* Hash alpha, beta from linear relation. */
    for (size_t i = 0; i < 2; i++) {
        blake3_hasher_update(&hasher, (const uint8_t *) alpha[i].data(),
                             8 * NTRU_DEGREE);
    }

    blake3_hasher_update(&hasher, (const uint8_t *) x.c1.data(), 8 * NTRU_DEGREE);
    blake3_hasher_update(&hasher, (const uint8_t *) y.c1.data(), 8 * NTRU_DEGREE);
    blake3_hasher_update(&hasher, (const uint8_t *) x.c2[0].data(),
                         8 * NTRU_DEGREE);
    blake3_hasher_update(&hasher, (const uint8_t *) y.c2[0].data(),
                         8 * NTRU_DEGREE);

    blake3_hasher_update(&hasher, (const uint8_t *) u.data(), 8 * NTRU_DEGREE);
    blake3_hasher_update(&hasher, (const uint8_t *) t.data(), 8 * NTRU_DEGREE);
    blake3_hasher_update(&hasher, (const uint8_t *) _t.data(), 8 * NTRU_DEGREE);

    blake3_hasher_finalize(&hasher, hash, BLAKE3_OUT_LEN);

    /* Sample challenge from RNG seeded with hash. */
    nfl::fastrandombytes_seed(hash);
    bdlop_sample_chal(beta);
    nfl::fastrandombytes_reseed();
}

static void poly_inverse(params::poly_q &inv, params::poly_q p) {
    std::array<mpz_t, params::poly_q::degree> coeffs;
    fmpz_t q;
    fmpz_mod_poly_t poly, irred;
    fmpz_mod_ctx_t ctx_q;

    fmpz_init(q);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
    }

    fmpz_set_mpz(q, params::poly_q::moduli_product());
    fmpz_mod_ctx_init(ctx_q, q);
    fmpz_mod_poly_init(poly, ctx_q);
    fmpz_mod_poly_init(irred, ctx_q);

    p.poly2mpz(coeffs);
    fmpz_mod_poly_set_coeff_ui(irred, params::poly_q::degree, 1, ctx_q);
    fmpz_mod_poly_set_coeff_ui(irred, 0, 1, ctx_q);

    for (size_t i = 0; i < params::poly_q::degree; i++) {
        flint_poly_set_coeff_mpz(poly, i, coeffs[i], ctx_q);
    }
    fmpz_mod_poly_invmod(poly, poly, irred, ctx_q);

    for (size_t i = 0; i < params::poly_q::degree; i++) {
        flint_poly_get_coeff_mpz(coeffs[i], poly, i, ctx_q);
    }

    inv.mpz2poly(coeffs);

    fmpz_clear(q);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
}

static void simul_inverse(params::poly_q inv[MSGS], params::poly_q m[MSGS]) {
    params::poly_q u;
    inv[0] = m[0];
    inv_tmp[0] = m[0];

    for (size_t i = 1; i < MSGS; i++) {
        inv_tmp[i] = m[i];
        inv[i] = inv[i - 1] * m[i];
    }

    u = inv[MSGS - 1];
    u.invntt_pow_invphi();
    poly_inverse(u, u);
    u.ntt_pow_phi();

    for (size_t i = MSGS - 1; i > 0; i--) {
        inv[i] = u * inv[i - 1];
        u = u * inv_tmp[i];
    }
    inv[0] = u;
}

static int rej_sampling(params::poly_q z[NTRU_WIDTH], params::poly_q v[NTRU_WIDTH],
                        uint64_t s2) {
    array<mpz_t, params::poly_q::degree> coeffs0, coeffs1;
    params::poly_q t;
    mpz_t dot, norm, qDivBy2, tmp;
    double r, M = 1.75;
    int64_t seed;
    mpf_t u;
    uint8_t buf[8];
    gmp_randstate_t state;
    int result;

    /// Constructors
    mpf_init(u);
    gmp_randinit_mt(state);
    mpz_inits(dot, norm, qDivBy2, tmp, nullptr);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs0[i], (params::poly_q::bits_in_moduli_product() << 2));
        mpz_init2(coeffs1[i], (params::poly_q::bits_in_moduli_product() << 2));
    }

    mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);
    mpz_set_ui(norm, 0);
    mpz_set_ui(dot, 0);
    for (int i = 0; i < NTRU_WIDTH; i++) {
        t = z[i];
        t.invntt_pow_invphi();
        t.poly2mpz(coeffs0);
        t = v[i];
        t.invntt_pow_invphi();
        t.poly2mpz(coeffs1);
        for (size_t i = 0; i < params::poly_q::degree; i++) {
            util::center(coeffs0[i], coeffs0[i],
                         params::poly_q::moduli_product(), qDivBy2);
            util::center(coeffs1[i], coeffs1[i],
                         params::poly_q::moduli_product(), qDivBy2);
            mpz_mul(tmp, coeffs0[i], coeffs1[i]);
            mpz_add(dot, dot, tmp);
            mpz_mul(tmp, coeffs1[i], coeffs1[i]);
            mpz_add(norm, norm, tmp);
        }
    }

    getrandom(buf, sizeof(buf), 0);
    memcpy(&seed, buf, sizeof(buf));
    gmp_randseed_ui(state, seed);
    mpf_urandomb(u, state, mpf_get_default_prec());

    r = -2.0 * mpz_get_d(dot) + mpz_get_d(norm);
    r = r / (2.0 * s2);
    r = exp(r) / M;
    result = mpf_get_d(u) > r;

    mpf_clear(u);
    mpz_clears(dot, norm, qDivBy2, tmp, nullptr);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs0[i]);
        mpz_clear(coeffs1[i]);
    }
    return result;
}

static void lin_prover(params::poly_q y[NTRU_WIDTH], params::poly_q _y[NTRU_WIDTH],
                       params::poly_q &t, params::poly_q &_t, params::poly_q &u,
                       commit_t x, commit_t _x, params::poly_q alpha[2],
                       comkey_t &key, vector<params::poly_q> r,
                       vector<params::poly_q> _r) {
    params::poly_q beta, tmp[NTRU_WIDTH], _tmp[NTRU_WIDTH];
    array<mpz_t, params::poly_q::degree> coeffs;
    mpz_t qDivBy2;
    int rej0, rej1;

    mpz_init(qDivBy2);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
    }
    mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);

    do {
        /* Prover samples y,y' from Gaussian. */
        for (int i = 0; i < NTRU_WIDTH; i++) {
            for (size_t k = 0; k < params::poly_q::degree; k++) {
                int64_t coeff = sample_z(0.0, NTRU_SIGMA_C);
                mpz_set_si(coeffs[k], coeff);
            }
            y[i].mpz2poly(coeffs);
            y[i].ntt_pow_phi();
            for (size_t k = 0; k < params::poly_q::degree; k++) {
                int64_t coeff = sample_z(0.0, NTRU_SIGMA_C);
                mpz_set_si(coeffs[k], coeff);
            }
            _y[i].mpz2poly(coeffs);
            _y[i].ntt_pow_phi();
        }

        t = y[0];
        _t = _y[0];
        for (int i = 0; i < NTRU_HEIGHT; i++) {
            for (int j = 0; j < NTRU_WIDTH - NTRU_HEIGHT; j++) {
                t = t + key.A1[i][j] * y[j + NTRU_HEIGHT];
                _t = _t + key.A1[i][j] * _y[j + NTRU_HEIGHT];
            }
        }

        u = 0;
        for (int i = 0; i < NTRU_WIDTH; i++) {
            u = u + alpha[0] * (key.A2[0][i] * y[i]) - (key.A2[0][i] * _y[i]);
        }

        /* Sample challenge. */
        lin_hash(beta, key, x, _x, alpha, u, t, _t);

        /* Prover */
        for (int i = 0; i < NTRU_WIDTH; i++) {
            tmp[i] = beta * r[i];
            _tmp[i] = beta * _r[i];
            y[i] = y[i] + tmp[i];
            _y[i] = _y[i] + _tmp[i];
        }
        rej0 = rej_sampling(y, tmp, NTRU_SIGMA_C * NTRU_SIGMA_C);
        rej1 = rej_sampling(_y, _tmp, NTRU_SIGMA_C * NTRU_SIGMA_C);
    } while (rej0 || rej1);

    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
    mpz_clear(qDivBy2);
}

static int lin_verifier(params::poly_q z[NTRU_WIDTH], params::poly_q _z[NTRU_WIDTH],
                        params::poly_q t, params::poly_q _t, params::poly_q u,
                        commit_t x, commit_t _x, params::poly_q alpha[2], comkey_t &key) {
    params::poly_q beta, v, _v, tmp;
    int result = 1;

    /* Sample challenge. */
    lin_hash(beta, key, x, _x, alpha, u, t, _t);

    /* Verifier checks norm, reconstruct from NTT representation. */
    for (int i = 0; i < NTRU_WIDTH; i++) {
        v = z[i];
        v.invntt_pow_invphi();
        result &= bdlop_test_norm(v, NTRU_SIGMA_C * NTRU_SIGMA_C);
        v = _z[i];
        v.invntt_pow_invphi();
        result &= bdlop_test_norm(v, NTRU_SIGMA_C * NTRU_SIGMA_C);
    }

    /* Verifier computes A1z and A1z'. */
    v = z[0];
    _v = _z[0];
    for (int i = 0; i < NTRU_HEIGHT; i++) {
        for (int j = 0; j < NTRU_WIDTH - NTRU_HEIGHT; j++) {
            v = v + key.A1[i][j] * z[j + NTRU_HEIGHT];
            _v = _v + key.A1[i][j] * _z[j + NTRU_HEIGHT];
        }
    }

    /* Being zero is preserved by the inverse transform, so these compare in
     * the NTT domain and skip it. */
    tmp = t + beta * x.c1 - v;
    result &= util::is_zero(tmp);
    tmp = _t + beta * _x.c1 - _v;
    result &= util::is_zero(tmp);

    v = 0;
    for (int i = 0; i < NTRU_WIDTH; i++) {
        v = v + alpha[0] * (key.A2[0][i] * z[i]) - (key.A2[0][i] * _z[i]);
    }
    t = (alpha[0] * x.c2[0] + alpha[1] - _x.c2[0]) * beta + u;

    result &= util::equal(t, v);
    return result;
}

void shuffle_hash(params::poly_q &beta, commit_t c[MSGS], commit_t d[MSGS],
                  params::poly_q _ms[MSGS], params::poly_q rho) {
    uint8_t hash[BLAKE3_OUT_LEN];
    blake3_hasher hasher;
    blake3_hasher_init(&hasher);

    blake3_hasher_init(&hasher);

    for (int i = 0; i < MSGS; i++) {
        blake3_hasher_update(&hasher, (const uint8_t *) _ms[i].data(),
                             8 * NTRU_DEGREE);
        blake3_hasher_update(&hasher, (const uint8_t *) c[i].c2[0].data(),
                             8 * NTRU_DEGREE);
        blake3_hasher_update(&hasher, (const uint8_t *) d[i].c2[0].data(),
                             8 * NTRU_DEGREE);
    }

    blake3_hasher_update(&hasher, (const uint8_t *) rho.data(),
                         8 * NTRU_DEGREE);
    blake3_hasher_finalize(&hasher, hash, BLAKE3_OUT_LEN);

    /* Sample challenge from RNG seeded with hash. */
    nfl::fastrandombytes_seed(hash);
    beta = nfl::uniform();
    nfl::fastrandombytes_reseed();
}

static void shuffle_prover(params::poly_q y[MSGS][NTRU_WIDTH],
                           params::poly_q _y[MSGS][NTRU_WIDTH], params::poly_q t[MSGS],
                           params::poly_q _t[MSGS], params::poly_q u[MSGS], commit_t d[MSGS],
                           params::poly_q s[MSGS], commit_t c[MSGS], params::poly_q ms[MSGS],
                           params::poly_q _ms[MSGS], vector<params::poly_q> r[MSGS],
                           params::poly_q rho, comkey_t &key) {
    params::poly_q t0;
    params::poly_q alpha[2], beta;

    /* Prover samples theta_i and computes commitments D_i. */
    for (size_t i = 0; i < MSGS - 1; i++) {
        theta[i] = nfl::ZO_dist();
        theta[i].ntt_pow_phi();
        if (i == 0) {
            t0 = theta[0] * _ms[0];
        } else {
            t0 = theta[i - 1] * ms[i] + theta[i] * _ms[i];
        }
        t0.invntt_pow_invphi();
        _r[i].resize(NTRU_WIDTH);
        bdlop_sample_rand(_r[i]);
        bdlop_commit(d[i], {t0}, key, _r[i]);
    }
    t0 = theta[MSGS - 2] * ms[MSGS - 1];
    t0.invntt_pow_invphi();
    _r[MSGS - 1].resize(NTRU_WIDTH);
    bdlop_sample_rand(_r[MSGS - 1]);
    bdlop_commit(d[MSGS - 1], {t0}, key, _r[MSGS - 1]);

    shuffle_hash(beta, c, d, _ms, rho);

    //Check relationship here
    simul_inverse(inv, _ms);
    for (size_t i = 0; i < MSGS - 1; i++) {
        if (i == 0) {
            s[0] = theta[0] * _ms[0] - beta * ms[0];
        } else {
            s[i] = theta[i - 1] * ms[i] + theta[i] * _ms[i] - s[i - 1] * ms[i];
        }
        s[i] = s[i] * inv[i];
    }

    /* Now run \Prod_LIN instances, one for each commitment. */
    for (size_t l = 0; l < MSGS; l++) {
        if (l < MSGS - 1) {
            t0 = s[l] * _ms[l];
        } else {
            if (MSGS & 1) {
                params::poly_q zero = 0;
                t0 = zero - beta * _ms[l];
            } else {
                t0 = beta * _ms[l];
            }
        }

        if (l == 0) {
            alpha[0] = beta;
        } else {
            alpha[0] = s[l - 1];
        }
        alpha[1] = t0;
        lin_prover(y[l], _y[l], t[l], _t[l], u[l], c[l], d[l], alpha, key, r[l],
                   _r[l]);
    }
}

static int shuffle_verifier(params::poly_q y[MSGS][NTRU_WIDTH],
                            params::poly_q _y[MSGS][NTRU_WIDTH], params::poly_q t[MSGS],
                            params::poly_q _t[MSGS], params::poly_q u[MSGS], commit_t d[MSGS],
                            params::poly_q s[MSGS], commit_t c[MSGS], params::poly_q _ms[MSGS],
                            params::poly_q rho, comkey_t &key) {
    params::poly_q alpha[2], beta;
    vector<params::poly_q> t0(1);
    int result = 1;

    shuffle_hash(beta, c, d, _ms, rho);
    for (size_t l = 0; l < MSGS; l++) {
        if (l < MSGS - 1) {
            t0[0] = s[l] * _ms[l];
        } else {
            if (MSGS & 1) {
                params::poly_q zero = 0;
                t0[0] = zero - beta * _ms[l];
            } else {
                t0[0] = beta * _ms[l];
            }
        }

        if (l == 0) {
            alpha[0] = beta;
        } else {
            alpha[0] = s[l - 1];
        }
        alpha[1] = t0[0];
        result &=
                lin_verifier(y[l], _y[l], t[l], _t[l], u[l], c[l], d[l], alpha,
                             key);
    }

    return result;
}

/**
 * Executes the run function of a cryptographic protocol involving commitments.
 *
 * @param com Array of commitments of size MSGS.
 * @param m Vector of vectors containing polynomial values.
 * @param _m Vector of vectors containing adjusted polynomial values.
 * @param key Commitment key.
 * @param r Array of vectors containing polynomial values of size MSGS.
 *
 * @return Returns the result of the shuffle_verifier function.
 *
 * @note This function extends commitments, adjusts keys, and then invokes
 *       the shuffle_prover and shuffle_verifier functions.
 */
static int run(commit_t com[MSGS], vector<params::poly_q> m,
               vector<params::poly_q> _m, comkey_t &key,
               vector<params::poly_q> r[MSGS]) {
    // Declare local variables for the function
    params::poly_q one, t1, rho;
    comkey_t _key;

    rho = nfl::uniform();  // Assign random values to rho
    for (size_t i = 0; i < MSGS; i++) {
        ms[i] = m[i];
        ms[i].ntt_pow_phi();  // Convert ms[i] to NTT domain
        ms[i] = ms[i] - rho;
        _ms[i] = _m[i];
        _ms[i].ntt_pow_phi();  // Convert _ms[i] to NTT domain
        _ms[i] = _ms[i] - rho;
        cs[i].c1 = com[i].c1;
        cs[i].c2 = com[i].c2;
        cs[i].c2[0] = cs[i].c2[0] - rho;
        t1 = m[i];
        t1.ntt_pow_phi();  // Convert t1 to NTT domain
        t1 = _m[i];
        t1.ntt_pow_phi();  // Convert t1 to NTT domain
    }

    // Call the shuffle_prover function
    shuffle_prover(y, _y, t, _t, u, d, s, cs, ms, _ms, r, rho, key);
    //return 1;

    // Call and return the result of the shuffle_verifier function
    return shuffle_verifier(y, _y, t, _t, u, d, s, cs, _ms, rho, key);
}

static void run2(commit_t com[MSGS], vector<params::poly_q> m,
               vector<params::poly_q> _m, comkey_t &key,
               vector<params::poly_q> r[MSGS]) {
    // Declare local variables for the function
    params::poly_q one, t1, rho;
    comkey_t _key;

    rho = nfl::uniform();  // Assign random values to rho
    for (size_t i = 0; i < MSGS; i++) {
        ms[i] = m[i];
        ms[i].ntt_pow_phi();  // Convert ms[i] to NTT domain
        ms[i] = ms[i] - rho;
        _ms[i] = _m[i];
        _ms[i].ntt_pow_phi();  // Convert _ms[i] to NTT domain
        _ms[i] = _ms[i] - rho;
        cs[i].c1 = com[i].c1;
        cs[i].c2 = com[i].c2;
        cs[i].c2[0] = cs[i].c2[0] - rho;
        t1 = m[i];
        t1.ntt_pow_phi();  // Convert t1 to NTT domain
        t1 = _m[i];
        t1.ntt_pow_phi();  // Convert t1 to NTT domain
    }

    // Call the shuffle_prover function
    shuffle_prover(y, _y, t, _t, u, d, s, cs, ms, _ms, r, rho, key);
}



#ifdef MAIN

static void test() {
    comkey_t key;
    vector<params::poly_q> m(MSGS), _m(MSGS);

    /* Generate commitment key-> */
    bdlop_keygen(key);
    for (int i = 0; i < MSGS; i++) {
        m[i] = nfl::ZO_dist();
        r[i].resize(NTRU_WIDTH);
        bdlop_sample_rand(r[i]);
        bdlop_commit(com[i], {m[i]}, key, r[i]);
    }

    /* Prover shuffles messages (only a circular shift for simplicity). */
    for (int i = 0; i < MSGS; i++) {
        _m[i] = m[(i + 1) % MSGS];
    }

    TEST_ONCE("polynomial inverse is correct")
    {
        params::poly_q alpha[2] = {nfl::uniform(), nfl::uniform()};

        poly_inverse(alpha[1], alpha[0]);
        alpha[0].ntt_pow_phi();
        alpha[1].ntt_pow_phi();
        alpha[0] = alpha[0] * alpha[1];
        alpha[0] = alpha[0] * alpha[1];
        alpha[0].invntt_pow_invphi();
        alpha[1].invntt_pow_invphi();
        TEST_ASSERT(util::equal(alpha[0], alpha[1]), end);
    }
    TEST_END;

    TEST_ONCE("shuffle proof is consistent")
    {
        TEST_ASSERT(run(com, m, _m, key, r) == 1, end);
    }
    TEST_END;

    end:
    return;
}

static void microbench() {
    params::poly_q alpha[2] = {nfl::uniform(), nfl::uniform()};

    alpha[0].ntt_pow_phi();
    alpha[1].ntt_pow_phi();

    BENCH_BEGIN("Polynomial addition")
        {
            BENCH_ADD(alpha[0] = alpha[0] + alpha[1]);
        }
    BENCH_END;

    BENCH_BEGIN("Polynomial multiplication")
        {
            BENCH_ADD(alpha[0] = alpha[0] * alpha[1]);
        }
    BENCH_END;

    alpha[0].invntt_pow_invphi();
    BENCH_BEGIN("Polynomial inverse")
        {
            BENCH_ADD(poly_inverse(alpha[1], alpha[0]));
        }
    BENCH_END;
}

static void bench() {
    comkey_t key;
    vector<params::poly_q> m(MSGS), _m(MSGS);
    params::poly_q y[NTRU_WIDTH], _y[NTRU_WIDTH], t, _t, u, alpha[2], beta;

    /* Generate commitment key-> */
    bdlop_keygen(key);
    for (int i = 0; i < MSGS; i++) {
        m[i] = nfl::ZO_dist();
        r[i].resize(NTRU_WIDTH);
        bdlop_sample_rand(r[i]);
        bdlop_commit(com[i], {m[i]}, key, r[i]);
    }

    /* Prover shuffles messages (only a circular shift for simplicity). */
    for (int i = 0; i < MSGS; i++) {
        _m[i] = m[(i + 1) % MSGS];
    }

    alpha[0] = nfl::ZO_dist();
    alpha[1] = nfl::ZO_dist();
    alpha[0].ntt_pow_phi();
    alpha[1].ntt_pow_phi();
    bdlop_sample_chal(beta);
    BENCH_BEGIN("linear hash")
        {
            BENCH_ADD(lin_hash(beta, key, com[0], com[1], alpha, u, t, _t));
        }
    BENCH_END;

    BENCH_BEGIN("linear proof")
        {
            BENCH_ADD(lin_prover(y, _y, t, _t, u, com[0], com[1], alpha, key, r[0],
                                 r[0]));
        }
    BENCH_END;

    BENCH_BEGIN("linear verifier")
        {
            BENCH_ADD(lin_verifier(y, _y, t, _t, u, com[0], com[1], alpha, key));
        }
    BENCH_END;

    BENCH_SMALL("shuffle-proof (N messages)", run(com, m, _m, key, r));
    BENCH_SMALL("shuffle proof (only proof)", run2(com, m, _m, key, r));
}

int main(int argc, char *argv[]) {
    shuffle_alloc();

    printf("\n** Tests for lattice-based shuffle proof:\n\n");
    test();

    printf("\n** Microbenchmarks for polynomial arithmetic:\n\n");
    microbench();

    printf("\n** Benchmarks for lattice-based shuffle proof:\n\n");
    bench();

    shuffle_free();
}

#endif