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


#define NTESTS 1

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
    PQCLEAN_MLKEM512_CLEAN_polyvec_frombytes(pk, packedpk);
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

static void make_rot_array(int16_t *p_o, poly *p_i){
    memcpy(p_o + KYBER_N, p_i, KYBER_N * sizeof(int16_t));

    for (unsigned int i = 0; i < KYBER_N; i++) {
        p_o[i] = -(p_i->coeffs[i]);
    }
}

static int pk_mask_check(uint8_t ct[KYBER_INDCPA_BYTES], const uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES]){
    unsigned int i, row, pos, pivot;
    uint8_t seed[KYBER_SYMBYTES];
    int16_t a0_inv;
    int16_t rot_ct[KYBER_N*2];
    polyvec pkpv, A_[KYBER_K], u;

    row = 0;
    pos = 0;
    pivot = 0;

    unpack_pk(&pkpv, seed, pk);
    PQCLEAN_MLKEM512_CLEAN_polyvec_decompress(&u, ct);               // Decompress ciphertext u
    PQCLEAN_MLKEM512_CLEAN_gen_matrix(A_, seed, 0);                 // at is NTT domain. No transpose
    PQCLEAN_MLKEM512_CLEAN_poly_invntt_tomont(&(A_[row].vec[pos])); // Convert A_ to INTT and Montgomery domain
    // TODO: 値を正にしとかないとだめ？
    // TODO: base処理
    a0_inv = 1306; // Inverse of 0x301 * 2^16 mod q

    for(int rot = 0; rot < KYBER_N; rot++) {
        make_rot_array(rot_ct, &u.vec[pos]);                            // Values range positive
    }

    print_poly((poly *)rot_ct);
    print_poly((poly *)(rot_ct + KYBER_N));


    return 0; // TODO: Implement the actual mask check logic
}

int main(void) {
    uint8_t key_a[CRYPTO_BYTES];
    // uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    // uint8_t ct[CRYPTO_CIPHERTEXTBYTES];
    // uint8_t sk_a[CRYPTO_SECRETKEYBYTES];
    // int i, j;

    // Read bytes from stdin
    uint8_t buffer[CRYPTO_CIPHERTEXTBYTES + CRYPTO_PUBLICKEYBYTES];
    size_t read_bytes = fread(buffer, 1, sizeof(buffer), stdin);

    // pk-mask check
    pk_mask_check(buffer, buffer + CRYPTO_CIPHERTEXTBYTES);

    // Decapsulation
    // crypto_kem_dec(key_a, buffer, buffer + CRYPTO_CIPHERTEXTBYTES);
    // printbytes(key_a, CRYPTO_BYTES);

    return 0;
}
