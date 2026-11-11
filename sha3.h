#ifndef SHA3_H
#define SHA3_H

#include <stdint.h>
#include <stddef.h>

typedef enum {
    SHA3_224 = 224,
    SHA3_256 = 256,
    SHA3_384 = 384,
    SHA3_512 = 512
} sha3_type;

typedef struct {
    uint64_t state[25];
    uint8_t  buffer[144];
    size_t   buffered;
    size_t   rate;
    size_t   digest_size;
} sha3_ctx;

int sha3_init(sha3_ctx *ctx, sha3_type type);

void sha3_update(
    sha3_ctx *ctx,
    const void *data,
    size_t len);

void sha3_finalize(
    sha3_ctx *ctx,
    uint8_t *digest);

int sha3(
    sha3_type type,
    const void *data,
    size_t len,
    uint8_t *digest);

size_t sha3_digest_size(sha3_type type);

#endif