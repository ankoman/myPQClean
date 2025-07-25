#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

//#include "api.h"
#include "randombytes.h"
#include "polyvec.h"
#include "poly.h"
#include "kem.h"
#include "indcpa.h"
#include "reduce.h"

#define N_coe 15

static void printbytes(const uint8_t *x, size_t xlen) {
    size_t i;
    for (i = 0; i < xlen; i++) {
        printf("%02x", x[i]);
    }
    printf("\n");
}

static void print_poly(const poly *p) {
    printf("poly.coeffs = [");
    for (size_t i = 0; i < KYBER_N; i++) {
        printf("%d", p->coeffs[i]);
        if (i < KYBER_N - 1) {
            printf(", ");
        }
    }
    printf("]\n");
}

static inline uint64_t rdtsc(){
    unsigned int lo, hi;
    __asm__ volatile (
        "cpuid\n"        // serialize
        "rdtsc\n"        // read time-stamp counter
        : "=a" (lo), "=d" (hi)
        : "a" (0)
        : "%ebx", "%ecx");
    return ((uint64_t)hi << 32) | lo;
}

// https://stackoverflow.com/a/1489985/1711232
#define PASTER(x, y) x##_##y
#define EVALUATOR(x, y) PASTER(x, y)
#define NAMESPACE(fun) EVALUATOR(PQCLEAN_NAMESPACE, fun)

#define CRYPTO_BYTES           NAMESPACE(CRYPTO_BYTES)
#define CRYPTO_PUBLICKEYBYTES  NAMESPACE(CRYPTO_PUBLICKEYBYTES)
#define CRYPTO_SECRETKEYBYTES  NAMESPACE(CRYPTO_SECRETKEYBYTES)
#define CRYPTO_CIPHERTEXTBYTES NAMESPACE(CRYPTO_CIPHERTEXTBYTES)

#define crypto_kem_keypair NAMESPACE(crypto_kem_keypair)
#define crypto_kem_keypair_derand NAMESPACE(crypto_kem_keypair_derand)
#define crypto_kem_enc NAMESPACE(crypto_kem_enc)
#define crypto_kem_dec NAMESPACE(crypto_kem_dec)

static void unpack_pk(polyvec *pk,
                      uint8_t seed[KYBER_SYMBYTES],
                      const uint8_t packedpk[KYBER_INDCPA_PUBLICKEYBYTES]) {
    PQCLEAN_MLKEM768_CLEAN_polyvec_frombytes(pk, packedpk);
    memcpy(seed, packedpk + KYBER_POLYVECBYTES, KYBER_SYMBYTES);
}

static int16_t comp(int16_t coeff) {
    uint64_t d0;
    int16_t t;
    t  = coeff;
    t += ((int16_t)t >> 15) & KYBER_Q;
    d0 = t;
    d0 <<= 10;
    d0 += 1665;
    d0 *= 1290167;
    d0 >>= 32;
    t = d0 & 0x3ff;
    return t;
}   

static int16_t decomp(uint16_t coeff) {
    return ((uint32_t)(coeff & 0x3FF) * KYBER_Q + 512) >> 10;
}   

static int16_t approx(int16_t coeff) {
    return decomp(comp(coeff));
}

static void make_rot_array_to_mont(int16_t *p_o, poly *p_i){
    // make all returned values positive
    int16_t t;
    int16_t r2 = 1353;
    
    for (unsigned int i = 0; i < KYBER_N; i++) {
        t = PQCLEAN_MLKEM768_CLEAN_montgomery_reduce((p_i->coeffs[i]) * r2);
        if (t < 0) {
            p_o[i + KYBER_N] = KYBER_Q + t;
            p_o[i] = -t;
        }else{
            p_o[i + KYBER_N] = t;
            p_o[i] = KYBER_Q - t;
        }
    }
}

static void make_rot_array(int16_t *p_o, int16_t *p_comp_o, poly *p_i){
    // make all returned values positive
    int16_t t;
    
    for (unsigned int i = 0; i < KYBER_N; i++) {
        t = p_i->coeffs[i];
        // Make negative
        // if (t > 1664){   
        //     t -= KYBER_Q;
        // }
        // p_o[i + KYBER_N] = t;
        // p_o[i] = -t;

        // Make positive
        if (t < 0) {
            p_o[i + KYBER_N] = KYBER_Q + t;
            p_o[i] = -t;
        }else{
            p_o[i + KYBER_N] = t;
            p_o[i] = KYBER_Q - t;
        }
        p_comp_o[i + KYBER_N] = comp(p_o[i + KYBER_N]);
        p_comp_o[i] = comp(p_o[i]);
    }
}
static void make_poly_positive(poly *p_i){
    // make all returned values positive
    int16_t t;
    
    for (unsigned int i = 0; i < KYBER_N; i++) {
        t = p_i->coeffs[i];
        // Make positive
        if (t < 0) {
            p_i->coeffs[i] = KYBER_Q + t;
        }
    }
}

static void polyvec_stay_comp(polyvec *r, const uint8_t a[KYBER_POLYVECCOMPRESSEDBYTES]) {
    unsigned int i, j, k;

    uint16_t t[4];
    for (i = 0; i < KYBER_K; i++) {
        for (j = 0; j < KYBER_N / 4; j++) {
            t[0] = (a[0] >> 0) | ((uint16_t)a[1] << 8);
            t[1] = (a[1] >> 2) | ((uint16_t)a[2] << 6);
            t[2] = (a[2] >> 4) | ((uint16_t)a[3] << 4);
            t[3] = (a[3] >> 6) | ((uint16_t)a[4] << 2);
            a += 5;

            for (k = 0; k < 4; k++) {
                r->vec[i].coeffs[4 * j + k] = t[k] & 0x3FF;
            }
        }
    }
}

static int cnt_invalid_coeffs(int16_t c, uint8_t base, int rot, poly *poly_pk, poly *poly_ct_comp){
    int cnt = 0;
    uint16_t val;

    int comp_base_idx = base + KYBER_N - rot;
    for (int i = 0; i < N_coe; i++) {
        val = comp((c * poly_pk->coeffs[base + i]) % KYBER_Q);
        if (val == poly_ct_comp->coeffs[comp_base_idx + i]) {
            cnt++;
        }
    }

    return cnt;
}

#define MULT t=(t*a)%KYBER_Q;
#define SQR  t=(t*t)%KYBER_Q;

static uint16_t inv(uint16_t a) {
    uint32_t t;
    t = a;
    SQR MULT SQR SQR SQR MULT SQR MULT SQR MULT
    SQR MULT SQR MULT SQR MULT SQR MULT SQR MULT
    return t;
}

static int pk_mask_check(uint8_t ct[KYBER_INDCPA_BYTES], const uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES]){
    unsigned int i, row, pos, base;
    uint8_t seed[KYBER_SYMBYTES], cnt;
    int16_t a0_inv, center, candidate, c;
    int16_t rot_ct[KYBER_N*2];
    int16_t rot_ct_comp[KYBER_N*2];
    polyvec pkpv, A_[KYBER_K], u;

    int DEBUG_ROT = 0, DEBUG_OFFSET = 0, DEBUG_CAND = 0, DEBUG_C = 0;

    row = 0;
    pos = 0;
    base = 0;

    unpack_pk(&pkpv, seed, pk);
    //polyvec_stay_comp(&u_comp, ct);
    PQCLEAN_MLKEM768_CLEAN_polyvec_decompress(&u, ct);              // Decompress ciphertext u
    PQCLEAN_MLKEM768_CLEAN_gen_matrix(A_, seed, 0);                 // at is NTT domain. No transpose
    PQCLEAN_MLKEM768_CLEAN_poly_invntt_tomont(&(A_[row].vec[pos])); // Convert A_ to INTT and Montgomery domain. Values can be negative
    make_poly_positive(&(A_[row].vec[pos]));                        // Make all values positive
    // TODO: 値を正にしとかないとだめ？
    // TODO: base処理
    // TODO: 0逆元の除外
    a0_inv = inv(A_[row].vec[pos].coeffs[0]); // Inverse of 0x301 * 2^16 mod q
    // printf("%d\n", a0_inv);
    
    make_rot_array(rot_ct, rot_ct_comp, &u.vec[pos]);                // u is decompressed. All returned values positive
    
    for(int rot = 0; rot < KYBER_N; rot++) {
        DEBUG_ROT++;
        center = rot_ct[base + KYBER_N - rot]; // Center value
        // TODO: Eliminate 0 center

        for (int offset = -2; offset <= 2; offset++) {
            DEBUG_OFFSET++;
            candidate = (center + offset) % KYBER_Q;         // TODO: Barrett reduction -> slow

            if (approx(candidate) == center){
                DEBUG_CAND++;
                c = (a0_inv * candidate) % KYBER_Q;
                // printf("c=%d\n", c);
                if(c < 416 || c > 2913){
                    DEBUG_C++;
                    cnt = cnt_invalid_coeffs(c, base, rot, &(A_[row].vec[pos]), (poly *)&rot_ct_comp);
                    // printf("cnt=%d\n", cnt);
                    if (cnt == N_coe) {
                        return 1; // pk-mask detected
                    }
                }
            }       
        }
    }
    // printf("%d, %d, %d, %d\n", DEBUG_ROT, DEBUG_OFFSET, DEBUG_CAND, DEBUG_C);
    // print_poly((poly *)&u_comp.vec[pos]);
    // print_poly((poly *)(rot_ct_comp));
    // print_poly((poly *)(rot_ct_comp + KYBER_N));


    return 0; // pk-mask undetected
}

#define N_TESTS 100000

#include <emmintrin.h>  // SSE2
#include <stdint.h>

static void flush_cache(void *p, size_t size) {
    char *cp = (char *)p;
    for (size_t i = 0; i < size; i += 64) {  // 64B = cache line size
        _mm_clflush(cp + i);
    }
}

int main(void) {
    uint8_t key_a[CRYPTO_BYTES];
    int res = 0;
    // uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    // uint8_t ct[CRYPTO_CIPHERTEXTBYTES];
    // uint8_t sk_a[CRYPTO_SECRETKEYBYTES];
    // int i, j;

    // Read bytes from stdin
    uint8_t buffer[CRYPTO_CIPHERTEXTBYTES + CRYPTO_PUBLICKEYBYTES];
    size_t read_bytes = fread(buffer, 1, sizeof(buffer), stdin);

    uint64_t start = rdtsc();
    // pk-mask check
    for (volatile int i = 0; i < N_TESTS; i++){
        //flush_cache(buffer, CRYPTO_CIPHERTEXTBYTES + CRYPTO_PUBLICKEYBYTES);
        res = res + pk_mask_check(buffer, buffer + CRYPTO_CIPHERTEXTBYTES);
        //crypto_kem_dec(key_a, buffer, buffer + CRYPTO_CIPHERTEXTBYTES);
    }
    uint64_t end = rdtsc();
    printf("Cycles: %f\n", (double)(end - start)/N_TESTS);
    printf("res= %d\n", res);

    // Decapsulation
    // crypto_kem_dec(key_a, buffer, buffer + CRYPTO_CIPHERTEXTBYTES);
    // printbytes(key_a, CRYPTO_BYTES);

    return res;
}
