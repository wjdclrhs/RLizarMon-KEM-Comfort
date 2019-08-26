#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "RLizarMon_KEM.h"


unsigned char pk[SEED_LEN+LWE_N];
unsigned char sk[LWE_N+LWE_N/8];

void Keygen_ring() {
	elapsed1 = clock();
	for (int l = 0; l < iter; ++l) {
		Keygen(pk, sk);
	}
	elapsed1 = clock() - elapsed1;

	printf("    Keygen Time: %f ms\n", elapsed1 * 1000. / CLOCKS_PER_SEC / iter);
}

void EncDecTest_RING() {
	// Set a messages
	unsigned char c[LWE_N+LWE_N+LWE_N/8] = {0,};
	unsigned char shared_k1[MESSAGE_LEN] = {0,};
	unsigned char shared_k2[MESSAGE_LEN] = {0,};
	int i, res = 0;

	elapsed1 = 0;
	elapsed2 = 0;

	for (int l = 0; l < iter; ++l) {
		for (i = 0; i < testnum; i++) {
			Enc(c, shared_k1, pk);
			res = Dec(shared_k2, c, sk, pk);
		}

		if (res == 1) {
			printf("    Decryption Validity Error Type 1 : c3 components\n");
			break;
		}

		if (res == 2) {
			printf("    Decryption Validity Error Type 2 : c2 components\n");
			break;
		}

		// Correctness check
		for (i = 0; i < MESSAGE_LEN; ++i) {
			if (shared_k1[i] != shared_k2[i]) {
				printf("Correctness Error, %d\n", i);
				break;
			}
		}
		if (i < MESSAGE_LEN) break;
	}
	printf("    Enc Time: %f ms\n", elapsed1 * 1000. / CLOCKS_PER_SEC / testnum / iter);
	printf("    Dec Time: %f ms\n", elapsed2 * 1000. / CLOCKS_PER_SEC / testnum / iter);






/*
	for (int l = 0; l < iter; ++l) {
		for (i = 0; i < testnum; i++) {
			Enc(c, m1, pk);
			res = Dec(m2, c, sk, pk);
			if (res == 1) {
				printf("    Decryption Validity Error Type 1 : c3 components\n");
				printf("iter: %d \t testnum: %d\n", iter, testnum);
				break;
			}

			if (res == 2) {
				printf("    Decryption Validity Error Type 2 : c2 components\n");
				printf("iter: %d \t testnum: %d\n", iter, testnum);
				break;
			}

			// Correctness check
			for (i = 0; i < MESSAGE_LEN; ++i) {
				if (m1[i] != m2[i]) {
					printf("Correctness Error, %d\n", i);
					break;
				}
			}
			if (i < MESSAGE_LEN) break;
		}
	}
	printf("    Enc Time: %f ms\n", elapsed1 * 1000. / CLOCKS_PER_SEC / testnum / iter);
	printf("    Dec Time: %f ms\n", elapsed2 * 1000. / CLOCKS_PER_SEC / testnum / iter);
*/
}


void main() {
	printf("\n  //////////////////////////////////////////////////////////////////\n\n");
	printf("\t\t"PARAMNAME" Parameter\n\n");
	printf("    LWE dimension: %d, \t\tLWR dimension: %d\n", LWE_N, LWE_N);
	printf("    Plaintext dimension: %d, \t\tPlaintext Modulus: %d bits\t\n", (MESSAGE_LEN*8), 1);
	printf("    Public Key modulus: %d bits, \tCiphertext modulus: %d bits\t\n\n", LOG_Q, LOG_P);
	printf("  //////////////////////////////////////////////////////////////////\n\n");
	printf("\t\t\tPerformance Test\n\n");

	// Key Generation
	Keygen_ring();

	// Enc and Dec
	EncDecTest_RING();
}
