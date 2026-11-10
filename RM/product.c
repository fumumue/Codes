
/*
 * concat2d.c
 *
 * Outer: 2D cyclic repetition code [4,1,4] over GF(2^20)
 * Inner: RM(1,4) [16,5,8] x Extended Hamming [8,4,4]
 * Product code: [128,20,32]
 *
 * Total binary code: [512,20,128]
 *
 * Inner decoder: exact nearest-neighbor decoding
 *
 * Compile:
 *   gcc -O3 -std=c11 -Wall -Wextra -o concat2d concat2d.c
 *
 * Run:
 *   ./concat2d
 *   ./concat2d 40
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define MSG_BITS 20
#define INNER_N 128
#define INNER_K 20
#define INNER_T 15
#define OUTER_N 4
#define TOTAL_N 512
#define TABLE_SIZE (1u << INNER_K)

typedef struct {
    uint64_t lo;
    uint64_t hi;
} Word128;

static Word128 *table;

/* XOR of two 128-bit words */
static Word128 wxor(Word128 a, Word128 b)
{
    Word128 r = {a.lo ^ b.lo, a.hi ^ b.hi};
    return r;
}

static int distance128(Word128 a, Word128 b)
{
    return __builtin_popcountll(a.lo ^ b.lo)
         + __builtin_popcountll(a.hi ^ b.hi);
}

static void setbit(Word128 *w, int pos, int value)
{
    if (!value) return;

    if (pos < 64)
        w->lo |= UINT64_C(1) << pos;
    else
        w->hi |= UINT64_C(1) << (pos - 64);
}

/* Extended Hamming [8,4,4] */
static void ham_encode(const int d[4], int c[8])
{
    memset(c, 0, 8 * sizeof(int));

    c[2] = d[0];
    c[4] = d[1];
    c[5] = d[2];
    c[6] = d[3];

    c[0] = d[0] ^ d[1] ^ d[3];
    c[1] = d[0] ^ d[2] ^ d[3];
    c[3] = d[1] ^ d[2] ^ d[3];

    for (int i = 0; i < 7; i++)
        c[7] ^= c[i];
}

/* RM(1,4) [16,5,8] */
static void rm_encode(const int a[5], int c[16])
{
    for (int x = 0; x < 16; x++) {
        int v = a[0];

        for (int j = 0; j < 4; j++)
            v ^= a[j + 1] & ((x >> j) & 1);

        c[x] = v;
    }
}

/*
 * Product encoder.
 *
 * 5 rows: RM information coefficients
 * 4 columns: Hamming information bits
 *
 * First encode Hamming horizontally.
 * Then encode RM vertically.
 */
static Word128 inner_encode(uint32_t msg)
{
    int matrix[5][8] = {{0}};
    Word128 out = {0, 0};

    for (int r = 0; r < 5; r++) {
        int d[4], h[8];

        for (int j = 0; j < 4; j++)
            d[j] = (msg >> (r * 4 + j)) & 1;

        ham_encode(d, h);

        for (int j = 0; j < 8; j++)
            matrix[r][j] = h[j];
    }

    for (int col = 0; col < 8; col++) {
        int a[5], rm[16];

        for (int r = 0; r < 5; r++)
            a[r] = matrix[r][col];

        rm_encode(a, rm);

        for (int r = 0; r < 16; r++)
            setbit(&out, r * 8 + col, rm[r]);
    }

    return out;
}

/* Build all 2^20 inner codewords */
static void build_table(void)
{
    table = malloc((size_t)TABLE_SIZE * sizeof(Word128));

    if (!table) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    Word128 basis[20];

    for (int i = 0; i < 20; i++)
        basis[i] = inner_encode(1u << i);

    table[0] = (Word128){0, 0};

    /*
     * Gray-independent recurrence:
     * code[m] = code[m without lowest set bit] XOR basis[bit]
     */
    for (uint32_t m = 1; m < TABLE_SIZE; m++) {
        unsigned bit = __builtin_ctz(m);
        table[m] = wxor(table[m & (m - 1)], basis[bit]);
    }

    printf("Inner codeword table: %u entries\n", TABLE_SIZE);
}

/*
 * Exact nearest-neighbor decoder.
 * Returns 1 if a unique codeword is within distance 15.
 */
static int inner_decode(Word128 received,
                        uint32_t *decoded,
                        int *best_distance)
{
    int best = INNER_N + 1;
    uint32_t best_msg = 0;

    for (uint32_t m = 0; m < TABLE_SIZE; m++) {
        int d = distance128(received, table[m]);

        if (d < best) {
            best = d;
            best_msg = m;

            if (best == 0)
                break;
        }
    }

    *best_distance = best;

    if (best > INNER_T)
        return 0;

    *decoded = best_msg;
    return 1;
}

/*
 * Outer 2D cyclic code:
 *
 *     a  a
 *     a  a
 *
 * Over GF(2^20).
 *
 * Generator polynomial:
 *     g(x,y) = 1 + x + y + xy
 *
 * Parameters: [4,1,4]
 */

/* Encode one symbol into four identical symbols */
static void outer_encode(uint32_t msg, uint32_t out[4])
{
    for (int i = 0; i < 4; i++)
        out[i] = msg;
}

/*
 * Outer bounded-distance decoder.
 *
 * erasure[i] = 1 means inner decoder failed.
 *
 * Unique decoding condition:
 *     2 * errors + erasures < 4
 */
static int outer_decode(const uint32_t in[4],
                        const int erasure[4],
                        uint32_t *decoded)
{
    int best_count = 0;
    uint32_t candidate = 0;
    int erasures = 0;

    for (int i = 0; i < 4; i++) {
        if (erasure[i]) {
            erasures++;
            continue;
        }

        int count = 0;

        for (int j = 0; j < 4; j++)
            if (!erasure[j] && in[i] == in[j])
                count++;

        if (count > best_count) {
            best_count = count;
            candidate = in[i];
        }
    }

    if (best_count == 0)
        return 0;

    int errors = 0;

    for (int i = 0; i < 4; i++)
        if (!erasure[i] && in[i] != candidate)
            errors++;

    if (2 * errors + erasures >= 4)
        return 0;

    *decoded = candidate;
    return 1;
}

int main(int argc, char **argv)
{
    int nerrors = 20;

    if (argc > 1)
        nerrors = atoi(argv[1]);

    if (nerrors < 0 || nerrors > TOTAL_N) {
        fprintf(stderr, "Errors must be 0..512\n");
        return 1;
    }

    build_table();

    srand(12345);

    uint32_t original = 0xABCDE;
    uint32_t outer[4];

    outer_encode(original, outer);

    Word128 received[4];

    for (int i = 0; i < 4; i++)
        received[i] = table[outer[i]];

    /* Generate distinct random error positions */
    int positions[TOTAL_N];

    for (int i = 0; i < TOTAL_N; i++)
        positions[i] = i;

    for (int i = TOTAL_N - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int t = positions[i];
        positions[i] = positions[j];
        positions[j] = t;
    }

    int block_errors[4] = {0};

    for (int i = 0; i < nerrors; i++) {
        int p = positions[i];
        int block = p / INNER_N;
        int bit = p % INNER_N;

        setbit(&received[block], bit, 0);

        /* Toggle the bit explicitly */
        if (bit < 64)
            received[block].lo ^= UINT64_C(1) << bit;
        else
            received[block].hi ^= UINT64_C(1) << (bit - 64);

        block_errors[block]++;
    }

    printf("\n=== Transmission ===\n");
    printf("Original message: %05X\n", original);
    printf("Injected errors: %d\n", nerrors);

    uint32_t inner_result[4] = {0};
    int erasure[4] = {0};

    printf("\n=== Inner decoding ===\n");

    for (int i = 0; i < 4; i++) {
        int dist = 0;

        int ok = inner_decode(received[i],
                              &inner_result[i],
                              &dist);

        erasure[i] = !ok;

        printf("Block %d: errors=%2d, nearest distance=%2d, ",
               i, block_errors[i], dist);

        if (ok)
            printf("decoded=%05X\n", inner_result[i]);
        else
            printf("ERASURE\n");
    }

    printf("\n=== Outer 2D cyclic decoding ===\n");

    uint32_t result = 0;

    int ok = outer_decode(inner_result, erasure, &result);

    if (ok) {
        printf("Decoded message: %05X\n", result);
        printf("Original message: %05X\n", original);

        if (result == original)
            printf("\nPASS: Correct recovery!\n");
        else
            printf("\nFAIL: Incorrect recovery!\n");
    } else {
        printf("Decoding failed.\n");
        printf("Too many inner errors/erasures.\n");
    }

    free(table);
    return 0;
}
