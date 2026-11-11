#include <stdint.h>
#include <stddef.h>
#include <string.h>

/* SHA-3 variants */
typedef enum {
    SHA3_224 = 224,
    SHA3_256 = 256,
    SHA3_384 = 384,
    SHA3_512 = 512
} sha3_type;

typedef struct {
    uint64_t state[25];
    uint8_t  buffer[144];  /* Maximum SHA-3 rate */
    size_t   buffered;
    size_t   rate;
    size_t   digest_size;
} sha3_ctx;

static const uint64_t keccak_rc[24] = {
    UINT64_C(0x0000000000000001),
    UINT64_C(0x0000000000008082),
    UINT64_C(0x800000000000808A),
    UINT64_C(0x8000000080008000),
    UINT64_C(0x000000000000808B),
    UINT64_C(0x0000000080000001),
    UINT64_C(0x8000000080008081),
    UINT64_C(0x8000000000008009),
    UINT64_C(0x000000000000008A),
    UINT64_C(0x0000000000000088),
    UINT64_C(0x0000000080008009),
    UINT64_C(0x000000008000000A),
    UINT64_C(0x000000008000808B),
    UINT64_C(0x800000000000008B),
    UINT64_C(0x8000000000008089),
    UINT64_C(0x8000000000008003),
    UINT64_C(0x8000000000008002),
    UINT64_C(0x8000000000000080),
    UINT64_C(0x000000000000800A),
    UINT64_C(0x800000008000000A),
    UINT64_C(0x8000000000008081),
    UINT64_C(0x8000000000008080),
    UINT64_C(0x0000000080000001),
    UINT64_C(0x8000000080008008)
};

/* Indexed by x + 5*y */
static const unsigned int keccak_rot[25] = {
     0,  1, 62, 28, 27,
    36, 44,  6, 55, 20,
     3, 10, 43, 25, 39,
    41, 45, 15, 21,  8,
    18,  2, 61, 56, 14
};

static uint64_t rol64(uint64_t x, unsigned int n)
{
    if (n == 0)
        return x;

    return (x << n) | (x >> (64 - n));
}

/* Keccak-f[1600], 24 rounds */
static void keccak_f1600(uint64_t a[25])
{
    uint64_t c[5], d[5], b[25];

    for (unsigned int round = 0; round < 24; ++round) {
        /* Theta */
        for (unsigned int x = 0; x < 5; ++x) {
            c[x] = a[x]
                 ^ a[x + 5]
                 ^ a[x + 10]
                 ^ a[x + 15]
                 ^ a[x + 20];
        }

        for (unsigned int x = 0; x < 5; ++x) {
            d[x] = c[(x + 4) % 5]
                 ^ rol64(c[(x + 1) % 5], 1);
        }

        for (unsigned int y = 0; y < 5; ++y) {
            for (unsigned int x = 0; x < 5; ++x) {
                a[x + 5 * y] ^= d[x];
            }
        }

        /* Rho and Pi */
        for (unsigned int y = 0; y < 5; ++y) {
            for (unsigned int x = 0; x < 5; ++x) {
                unsigned int nx = y;
                unsigned int ny = (2 * x + 3 * y) % 5;

                b[nx + 5 * ny] =
                    rol64(a[x + 5 * y], keccak_rot[x + 5 * y]);
            }
        }

        /* Chi */
        for (unsigned int y = 0; y < 5; ++y) {
            for (unsigned int x = 0; x < 5; ++x) {
                a[x + 5 * y] =
                    b[x + 5 * y]
                    ^ ((~b[(x + 1) % 5 + 5 * y])
                       & b[(x + 2) % 5 + 5 * y]);
            }
        }

        /* Iota */
        a[0] ^= keccak_rc[round];
    }
}

/* XOR one rate-sized block into the state, then permute */
static void sha3_absorb_block(
    sha3_ctx *ctx,
    const uint8_t *block)
{
    for (size_t i = 0; i < ctx->rate; ++i) {
        size_t lane = i / 8;
        unsigned int shift = (unsigned int)(8 * (i % 8));

        ctx->state[lane] ^= (uint64_t)block[i] << shift;
    }

    keccak_f1600(ctx->state);
}

/*
 * Initialize the context.
 * Returns 0 on success, -1 for an unsupported SHA-3 type.
 */
int sha3_init(sha3_ctx *ctx, sha3_type type)
{
    size_t rate;
    size_t digest_size;

    switch (type) {
    case SHA3_224:
        rate = 144;
        digest_size = 28;
        break;

    case SHA3_256:
        rate = 136;
        digest_size = 32;
        break;

    case SHA3_384:
        rate = 104;
        digest_size = 48;
        break;

    case SHA3_512:
        rate = 72;
        digest_size = 64;
        break;

    default:
        return -1;
    }

    memset(ctx, 0, sizeof(*ctx));
    ctx->rate = rate;
    ctx->digest_size = digest_size;

    return 0;
}

/* Add input bytes. Input may be split across multiple calls. */
void sha3_update(
    sha3_ctx *ctx,
    const void *data,
    size_t len)
{
    const uint8_t *input = (const uint8_t *)data;

    while (len > 0) {
        size_t available = ctx->rate - ctx->buffered;
        size_t take = (len < available) ? len : available;

        memcpy(ctx->buffer + ctx->buffered, input, take);

        ctx->buffered += take;
        input += take;
        len -= take;

        if (ctx->buffered == ctx->rate) {
            sha3_absorb_block(ctx, ctx->buffer);
            ctx->buffered = 0;
        }
    }
}

/*
 * Finish hashing and write the digest.
 * The output buffer must have at least the digest size for this variant.
 */
void sha3_finalize(
    sha3_ctx *ctx,
    uint8_t *digest)
{
    /* SHA-3 domain separation and multi-rate padding */
    memset(ctx->buffer + ctx->buffered,
           0,
           ctx->rate - ctx->buffered);

    ctx->buffer[ctx->buffered] = 0x06;
    ctx->buffer[ctx->rate - 1] |= 0x80;

    sha3_absorb_block(ctx, ctx->buffer);

    /* Extract bytes in little-endian lane order */
    for (size_t i = 0; i < ctx->digest_size; ++i) {
        size_t lane = i / 8;
        unsigned int shift = (unsigned int)(8 * (i % 8));

        digest[i] = (uint8_t)(ctx->state[lane] >> shift);
    }

    memset(ctx, 0, sizeof(*ctx));
}

/* One-shot interface */
int sha3(
    sha3_type type,
    const void *data,
    size_t len,
    uint8_t *digest)
{
    sha3_ctx ctx;

    if (sha3_init(&ctx, type) != 0)
        return -1;

    sha3_update(&ctx, data, len);
    sha3_finalize(&ctx, digest);

    return 0;
}

/* Return digest length in bytes, or 0 for an invalid type */
size_t sha3_digest_size(sha3_type type)
{
    switch (type) {
    case SHA3_224: return 28;
    case SHA3_256: return 32;
    case SHA3_384: return 48;
    case SHA3_512: return 64;
    default:       return 0;
    }
}
