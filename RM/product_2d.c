
/*
 * concatenated_2d_product.c
 *
 * Inner code : binary product code [9,4,4]
 * Outer code : 2D cyclic code over GF(16), [9,4,4]
 * Total      : binary concatenated code [81,16,>=16]
 *
 * gcc -O2 -std=c11 -o concat concatenated_2d_product.c
 * ./concat
 */

#include <stdio.h>
#include <stdint.h>
#include <string.h>

#define Q 16
#define INNER_N 9
#define INNER_K 4
#define OUTER_N 9
#define OUTER_K 4
#define TOTAL_N 81
#define TOTAL_K 16

/* Test errors: one bit in each of three inner blocks */
static const int ERROR_POSITIONS[] = {0, 18, 36};
#define ERROR_COUNT 3

/* =========================================
   GF(16), primitive polynomial x^4+x+1
   ========================================= */

static uint8_t gf_mul(uint8_t a, uint8_t b)
{
    uint8_t r = 0;

    while (b) {
        if (b & 1)
            r ^= a;

        b >>= 1;
        a <<= 1;

        if (a & 0x10)
            a ^= 0x13;
    }

    return r & 15;
}

/* =========================================
   INNER CODE: binary [9,4,4]
   ========================================= */

static uint16_t inner_encode(uint8_t value)
{
    int d[4];

    for (int i = 0; i < 4; i++)
        d[i] = (value >> (3-i)) & 1;

    int c[3][3] = {{0}};

    c[0][0] = d[0];
    c[0][1] = d[1];
    c[1][0] = d[2];
    c[1][1] = d[3];

    /* Row parity */
    c[0][2] = c[0][0] ^ c[0][1];
    c[1][2] = c[1][0] ^ c[1][1];

    /* Column parity */
    c[2][0] = c[0][0] ^ c[1][0];
    c[2][1] = c[0][1] ^ c[1][1];

    /* Corner parity */
    c[2][2] = c[0][2] ^ c[1][2];

    uint16_t result = 0;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            result = (result << 1) | c[i][j];

    return result;
}

static int popcount9(uint16_t x)
{
    int n = 0;

    while (x) {
        n += x & 1;
        x >>= 1;
    }

    return n;
}

static int inner_decode(
    uint16_t received,
    uint8_t *decoded,
    int *distance)
{
    int best = 10;
    int ties = 0;
    uint8_t answer = 0;

    for (int m = 0; m < 16; m++) {
        uint16_t candidate = inner_encode(m);

        int d = popcount9(received ^ candidate);

        if (d < best) {
            best = d;
            answer = (uint8_t)m;
            ties = 1;
        }
        else if (d == best) {
            ties++;
        }
    }

    *decoded = answer;
    *distance = best;

    return ties == 1;
}

/* =========================================
   OUTER CODE: 2D cyclic [9,4,4] over GF(16)
   
   Generator:
       g(x,y) = (1+x)(1+y)

   Ring:
       GF(16)[x,y]/(x^3-1,y^3-1)
   ========================================= */

static void outer_encode(
    const uint8_t message[4],
    uint8_t codeword[9])
{
    uint8_t u[3][3] = {{0}};
    uint8_t c[3][3] = {{0}};

    u[0][0] = message[0];
    u[0][1] = message[1];
    u[1][0] = message[2];
    u[1][1] = message[3];

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {

            int im = (i + 2) % 3;
            int jm = (j + 2) % 3;

            c[i][j] =
                u[i][j]
                ^ u[im][j]
                ^ u[i][jm]
                ^ u[im][jm];
        }
    }

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            codeword[3*i+j] = c[i][j];
}

/* =========================================
   OUTER DECODER
   Exhaustive nearest-codeword decoding
   16^4 = 65536 candidates
   ========================================= */

static int outer_decode(
    const uint8_t received[9],
    uint8_t decoded[4],
    int *distance)
{
    uint8_t best_message[4] = {0};
    uint8_t candidate[9];

    int best = 10;
    int ties = 0;

    for (int m = 0; m < 65536; m++) {

        uint8_t message[4] = {
            (m >> 12) & 15,
            (m >> 8) & 15,
            (m >> 4) & 15,
            m & 15
        };

        outer_encode(message, candidate);

        int d = 0;

        for (int i = 0; i < 9; i++)
            d += (received[i] != candidate[i]);

        if (d < best) {
            best = d;
            ties = 1;
            memcpy(best_message, message, 4);
        }
        else if (d == best) {
            ties++;
        }

        if (best == 0)
            break;
    }

    memcpy(decoded, best_message, 4);
    *distance = best;

    return ties == 1;
}

/* =========================================
   CONCATENATED ENCODER
   ========================================= */

static void encode(
    const uint8_t message[4],
    uint8_t transmitted[81])
{
    uint8_t outer[9];

    outer_encode(message, outer);

    for (int i = 0; i < 9; i++) {

        uint16_t inner = inner_encode(outer[i]);

        for (int j = 0; j < 9; j++) {
            transmitted[i*9+j] =
                (inner >> (8-j)) & 1;
        }
    }
}

/* =========================================
   CONCATENATED DECODER
   ========================================= */

static int decode(
    const uint8_t received[81],
    uint8_t decoded[4])
{
    uint8_t inner_symbols[9];

    printf("\nINNER DECODING\n");

    for (int i = 0; i < 9; i++) {

        uint16_t block = 0;

        for (int j = 0; j < 9; j++)
            block = (block << 1) |
                    received[i*9+j];

        int distance;
        int unique;

        unique = inner_decode(
            block,
            &inner_symbols[i],
            &distance
        );

        printf(
            "Block %d: symbol=%X distance=%d unique=%s\n",
            i,
            inner_symbols[i],
            distance,
            unique ? "YES" : "NO"
        );
    }

    int outer_distance;

    printf("\nOUTER DECODING\n");

    int unique = outer_decode(
        inner_symbols,
        decoded,
        &outer_distance
    );

    printf("Outer symbol errors: %d\n", outer_distance);
    printf("Unique decoding: %s\n",
           unique ? "YES" : "NO");

    return unique;
}

/* =========================================
   DISPLAY
   ========================================= */

static void print_bits(const char *name,
                       const uint8_t bits[81])
{
    printf("%s\n", name);

    for (int i = 0; i < 81; i++) {
        printf("%d", bits[i]);

        if ((i+1) % 9 == 0)
            printf(" ");
    }

    printf("\n");
}

/* =========================================
   MAIN TEST
   ========================================= */

int main(void)
{
    uint8_t message[4] = {
        0xB, 0x6, 0xC, 0x3
    };

    uint8_t transmitted[81];
    uint8_t received[81];
    uint8_t decoded[4];

    encode(message, transmitted);

    memcpy(received, transmitted, 81);

    /* Inject errors */

    for (int i = 0; i < ERROR_COUNT; i++) {
        int p = ERROR_POSITIONS[i];
        received[p] ^= 1;
    }

    printf("================================\n");
    printf("2D CYCLIC + PRODUCT CONCATENATED CODE\n");
    printf("================================\n");

    printf("Field: GF(16)\n");
    printf("Inner: [9,4,4] binary product code\n");
    printf("Outer: [9,4,4] 2D cyclic code\n");
    printf("Total: [81,16,>=16] binary code\n");

    printf("\nOriginal message: ");

    for (int i = 0; i < 4; i++)
        printf("%X ", message[i]);

    printf("\n");

    print_bits("Transmitted:", transmitted);
    print_bits("Received:", received);

    printf("\nERROR VECTOR\n");

    uint8_t error_vector[81];

    for (int i = 0; i < 81; i++)
        error_vector[i] = transmitted[i] ^ received[i];

    print_bits("Injected errors:", error_vector);

    printf("Error positions:");

    for (int i = 0; i < ERROR_COUNT; i++)
        printf(" %d", ERROR_POSITIONS[i]);

    printf("\nError weight: %d\n", ERROR_COUNT);

    if (!decode(received, decoded)) {
        printf("\nDECODING AMBIGUOUS\n");
        return 1;
    }

    printf("\n================================\n");
    printf("FINAL RESULT\n");
    printf("================================\n");

    printf("Original : ");

    for (int i = 0; i < 4; i++)
        printf("%X ", message[i]);

    printf("\nRecovered: ");

    for (int i = 0; i < 4; i++)
        printf("%X ", decoded[i]);

    printf("\n");

    uint8_t corrected[81];

    encode(decoded, corrected);

    int success = memcmp(
        transmitted, corrected, 81
    ) == 0;

    printf("SUCCESS: %s\n",
           success ? "YES" : "NO");

    printf("\nRecovered error vector:\n");

    for (int i = 0; i < 81; i++)
        printf("%d", received[i] ^ corrected[i]);

    printf("\n");

    return success ? 0 : 1;
}
