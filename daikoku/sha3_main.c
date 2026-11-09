#include "sha3.h"
#include <stdio.h>
#include <stdint.h>

int main(void)
{
    uint8_t digest[64];

    sha3(SHA3_512, "abc", 3, digest);

    for (size_t i = 0; i < sha3_digest_size(SHA3_512); ++i)
        printf("%02x", digest[i]);

    putchar('\n');
    return 0;
}

