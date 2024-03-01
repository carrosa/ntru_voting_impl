#include <cstddef>

#include <gmpxx.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <memory>
#include <nfl.hpp>
#include <cmath>

#include "bench.h"
#include <sys/random.h>

using namespace std;

#ifndef COMMON_H
#define COMMON_H

/* Parameter v in the commitment  scheme (laximum l1-norm of challs). */
#define NONZERO     36
/* Width k of the commitment matrix. */
#define WIDTH        4
/* Height of the commitment matrix. */
#define HEIGHT        1
/* Dimension of the committed messages. */
#ifndef SIZE
#define SIZE        1
#endif
/* Degree of the irreducible polynomial. */
#define DEGREE      2048//4096
/* Sigma for the commitment gaussian distribution. */
#define SIGMA_C     (1u << 12)
/* Sigma for the boundness proof. */
/* Parties that run the distributed decryption protocol. */
#define PARTIES     4
/* Security level for Distributed Decryption. */

/*
 * NTRU DEFINEs
 * */
#define NTRU_PRIMEQ 576460752303439873
#define NTRU_PRIMEP 2
#define NTRU_SIGMA 7.12
#define NTRU_DEGREE 2048
#define NTRU_BOUND_D "144183102358236620"
#define NTRU_PARTIES 4
/* Dimension of the committed messages. */
#ifndef NTRU_SIZE
#define NTRU_SIZE 1
#endif
/* Width k of the commitment matrix. */
#define NTRU_WIDTH 4
/* Height of the commitment matrix. */
#define NTRU_HEIGHT 1

#define NTRU_SIGMA_C (1u << 10)
#define NTRU_NONZERO 14


namespace params {
    using poly_p = nfl::poly_from_modulus<uint32_t, DEGREE, 30>;
    using poly_q = nfl::poly_from_modulus<uint64_t, DEGREE, 62>;
    using poly_big = nfl::poly_from_modulus<uint64_t, 4 * DEGREE, 62>;
}

namespace ntru_params {
    using poly_p = nfl::poly_from_modulus<uint32_t, NTRU_DEGREE, 30>;
    using poly_q = nfl::poly_from_modulus<uint64_t, NTRU_DEGREE, 62>;
    using poly_big = nfl::poly_from_modulus<uint64_t, 4 * DEGREE, 62>;
}

/*============================================================================*/
/* Type definitions                                                           */
/*============================================================================*/

/* Class that represents a commitment key pair. */
class comkey_t {
public:
    params::poly_q A1[HEIGHT][WIDTH - HEIGHT];
    params::poly_q A2[SIZE][WIDTH];
};

class ntru_comkey_t {
public:
    ntru_params::poly_q A1[HEIGHT][WIDTH - HEIGHT];
    ntru_params::poly_q A2[WIDTH];
};

/* Class that represents a commitment in CRT representation. */
class commit_t {
public:
    params::poly_q c1;
    vector<params::poly_q> c2;
};

class ntru_commit_t {
public:
    ntru_params::poly_q c1;
    ntru_params::poly_q c2;
};

/* Class that represents a BGV key pair. */
class bgvkey_t {
public:
    params::poly_q a;
    params::poly_q b;
};

class bgvenc_t {
public:
    params::poly_q u;
    params::poly_q v;
};

#include "util.hpp"

void bdlop_sample_rand(vector<params::poly_q> &r);

void bdlop_sample_chal(params::poly_q &f);

bool bdlop_test_norm(params::poly_q r, double_t sigma_sqr);

void bdlop_commit(commit_t &com, vector<params::poly_q> m, comkey_t &key, vector<params::poly_q> r);

int bdlop_open(commit_t &com, vector<params::poly_q> m, comkey_t &key, vector<params::poly_q> r, params::poly_q &f);

void bdlop_keygen(comkey_t &key);

void bgv_sample_message(params::poly_p &r);

void bgv_sample_short(params::poly_q &r);

void bgv_keygen(bgvkey_t &pk, params::poly_q &sk);

void bgv_encrypt(bgvenc_t &c, bgvkey_t &pk, params::poly_p &m);

void bgv_decrypt(params::poly_p &m, bgvenc_t &c, params::poly_q &sk);

// ntru
void ntru_bdlop_sample_rand(vector<ntru_params::poly_q> &r);

void ntru_bdlop_sample_chal(ntru_params::poly_q &f);

bool ntru_bdlop_test_norm(ntru_params::poly_q r, double_t sigma_sqr);

void ntru_bdlop_commit(ntru_commit_t &com, ntru_params::poly_q &m, ntru_comkey_t &key, vector<ntru_params::poly_q> r);

int ntru_bdlop_open(ntru_commit_t &com, ntru_params::poly_q m, ntru_comkey_t &key, vector<ntru_params::poly_q> r, ntru_params::poly_q &f);

void ntru_bdlop_keygen(ntru_comkey_t &key);

void ntru_keygen(ntru_params::poly_q &pk, ntru_params::poly_q &sk);
void ntru_sample_message(ntru_params::poly_p &r);
void ntru_encrypt(ntru_params::poly_q &c, ntru_params::poly_q &pk, ntru_params::poly_p &m);


#endif
