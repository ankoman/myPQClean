#ifndef PQCLEAN_MLDSA44_AARCH64_API_H
#define PQCLEAN_MLDSA44_AARCH64_API_H

/*
 * This file is dual licensed
 * under Apache 2.0 (https://www.apache.org/licenses/LICENSE-2.0.html)
 * or public domain at https://github.com/pq-crystals/dilithium
 */

#include <stddef.h>
#include <stdint.h>

#define PQCLEAN_MLDSA44_AARCH64_CRYPTO_PUBLICKEYBYTES 1312
#define PQCLEAN_MLDSA44_AARCH64_CRYPTO_SECRETKEYBYTES 2560
#define PQCLEAN_MLDSA44_AARCH64_CRYPTO_BYTES          2420
#define PQCLEAN_MLDSA44_AARCH64_CRYPTO_ALGNAME       "ML-DSA-44"

int PQCLEAN_MLDSA44_AARCH64_crypto_sign_keypair(uint8_t *pk, uint8_t *sk);

int PQCLEAN_MLDSA44_AARCH64_crypto_sign_signature_ctx(uint8_t *sig, size_t *siglen,
        const uint8_t *m, size_t mlen,
        const uint8_t *ctx, size_t ctxlen,
        const uint8_t *sk);

int PQCLEAN_MLDSA44_AARCH64_crypto_sign_ctx(uint8_t *sm, size_t *smlen,
        const uint8_t *m, size_t mlen,
        const uint8_t *ctx, size_t ctxlen,
        const uint8_t *sk);

int PQCLEAN_MLDSA44_AARCH64_crypto_sign_verify_ctx(const uint8_t *sig, size_t siglen,
        const uint8_t *m, size_t mlen,
        const uint8_t *ctx, size_t ctxlen,
        const uint8_t *pk);

int PQCLEAN_MLDSA44_AARCH64_crypto_sign_open_ctx(uint8_t *m, size_t *mlen,
        const uint8_t *sm, size_t smlen,
        const uint8_t *ctx, size_t ctxlen,
        const uint8_t *pk);

int PQCLEAN_MLDSA44_AARCH64_crypto_sign_signature(uint8_t *sig, size_t *siglen,
        const uint8_t *m, size_t mlen,
        const uint8_t *sk);

int PQCLEAN_MLDSA44_AARCH64_crypto_sign(uint8_t *sm, size_t *smlen,
                                        const uint8_t *m, size_t mlen,
                                        const uint8_t *sk);

int PQCLEAN_MLDSA44_AARCH64_crypto_sign_verify(const uint8_t *sig, size_t siglen,
        const uint8_t *m, size_t mlen,
        const uint8_t *pk);

int PQCLEAN_MLDSA44_AARCH64_crypto_sign_open(uint8_t *m, size_t *mlen,
        const uint8_t *sm, size_t smlen,
        const uint8_t *pk);

#endif
