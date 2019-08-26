#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "RLizarMon_KEM.h"
#include "randombytes.h"
#include "fips202.h"
#include "bch.h"
#include "ecc.h"

int Keygen(unsigned char *pk, unsigned char *sk){
	unsigned char pk_a[LWE_N];
	unsigned char pk_b[LWE_N*2];
	//unsigned char pk[SEED_LEN+LWE_N]; // seed_a + pk_b
	//unsigned char sk[LWE_N] = {0,};
	unsigned char seed_a[SEED_LEN];
	uint16_t sk_s[HS];
	int i, j;

    // Gen poly a
	randombytes(seed_a, SEED_LEN);	
	shake256(pk_a, LWE_N, seed_a, SEED_LEN);
/*	for(i = 0; i < LWE_N; ++i){
		printf("pk_a[%d]: %d \n", i, pk_a[i]);
	}
*/

    // gen poly s
	memset(sk, 0, LWE_N);
	unsigned char seed_s[HS*4];
	unsigned int sk_random_idx;
	int hw, count = 0;

	randombytes(seed_s, HS*4);

	while (hw < HS) {
	sk_random_idx = seed_s[count++]; 
	sk_random_idx <<= 8;
	sk_random_idx ^= seed_s[count++];
	sk_random_idx &= (LWE_N - 1); // (seed1) || (seed2) 
	if (sk[sk_random_idx] == 0) {
		sk[sk_random_idx] = (seed_s[count++] & 0x02) - 1;
		hw++;
	}
	if (count >= HS*4 - 3) {
		//printf("New seed_s \n");
		randombytes(seed_s, HS*4);
		count = 0;
	}
	}

//	int one=0;
//	int minus=0;
//	int zero=0;
//	for(i = 0; i<LWE_N; ++i){
//		printf("sk[%d]: %d \n", i, sk[i]);
//		if (sk_s[i] == 0x01){
//			one++;
//		}
//		else if (sk_s[i] == 0xff){
//			minus++;
//		}
//		else { zero++;}
//	}
//	printf("one: %d \t minus: %d \t zero: %d \n", one, minus, zero);



  // gen s_idx
	int neg_start = 0, back_position = HS;	
	for (i = 0; i < LWE_N; ++i) {
		if (sk[i] == 0x01){ 
			sk_s[neg_start++] = i;
		}
		if (sk[i] == 0xff){
			sk_s[--back_position] = i;
		}
	}
/*	
	for(i = 0; i<HS; ++i){
		printf("sk_s[%d]: %d \n", i, sk_s[i]);
	}
*/

	// gen u
	unsigned char u[LWE_N/8];
	randombytes(u, LWE_N/8);


	//concat sk = sk || u
	for (i=LWE_N; i<LWE_N+LWE_N/8; ++i){
		sk[i] = u[i-LWE_N];
//		printf("sk[%d]: %d \n", i, sk[i]);
	}





	// Initialize b as an error polynomial e
	unsigned char b0, b1, tmp2[LWE_N/4];
	memset(pk_b, 0, LWE_N*2);
	randombytes(tmp2,LWE_N/4); // tmp2[0]'s 0, 1bit = pk_b[0], tmp2[0]'s 2,3bit = pk_b[1]and so on
    // Centered Binomial Distribution
	for(j=0; j<LWE_N/4; ++j){
		b0 = tmp2[j] & 0x01;
		tmp2[j] = tmp2[j] >> 1;
		b1 = tmp2[j] & 0x01;
		pk_b[j*4+0] = b0 -b1;
		tmp2[j] = tmp2[j] >> 1;
		b0 = tmp2[j] & 0x01;
		tmp2[j] = tmp2[j] >> 1;
		b1 = tmp2[j] & 0x01;
		pk_b[j*4+1] = b0 -b1;
		tmp2[j] = tmp2[j] >> 1;
		b0 = tmp2[j] & 0x01;
		tmp2[j] = tmp2[j] >> 1;
		b1 = tmp2[j] & 0x01;
		pk_b[j*4+2] = b0 -b1;
		tmp2[j] = tmp2[j] >> 1;
		b0 = tmp2[j] & 0x01;
		tmp2[j] = tmp2[j] >> 1;
		b1 = tmp2[j] & 0x01;
		pk_b[j*4+3] = b0 -b1;
		tmp2[j] = tmp2[j] >> 1;		
	}

//	int one=0;
//	int minus=0;
//	int zero=0;
//	for(i = 0; i<LWE_N; ++i){
//		printf("pk_b[%d]: %d \n", i, pk_b[i]);
//		if (pk_b[i] == 0x01){
//			one++;
//		}
//		else if (pk_b[i] == 0xff){
//			minus++;
//		}
//		else { zero++;}
//	}
//	printf("one: %d \t minus: %d \t zero: %d \n", one, minus, zero);





    //mult a*s and add e(pk_b)
	for (i = 0; i < HS; ++i) {
		uint16_t deg = sk_s[i];
		uint16_t branch = (2 * ((i - neg_start) >> sft & 0x1) - 1);
		for (int j = 0; j < LWE_N; ++j) {
			pk_b[deg + j] -= branch * pk_a[j];  
		}
	}
	for (j = 0; j < LWE_N; ++j) {
		pk_b[j] -= pk_b[LWE_N + j];
	}
/*
	for (i=0;i<LWE_N;++i){
		printf("pk_b[%d]: %d \n", i, pk_b[i]);
	}
*/

    // Concat seed_genA || pk_b
	for (i = 0; i < SEED_LEN; ++i) {pk[i] = seed_a[i];}
	for (i = 0; i < LWE_N; ++i) {pk[SEED_LEN + i] = pk_b[i];}

	return 0;
}







int Enc(unsigned char *c, unsigned char *shared_k, const unsigned char *pk){
	start = clock();
	int i, j = 0;
	unsigned char c1h_a[LWE_N*2]={0,};
	unsigned char c1h_b[LWE_N*2]={0,};	
	unsigned char *hash = NULL;
	unsigned char *hash_t = NULL;

	//c[0]~[31] : c1, c[32]~[543]: c2, c[544]~[1055]: c3, c[1056]~[1119]: c4

	// Generate a random polynomial delta
	unsigned char delta[size_of_delta];
	randombytes(delta, (size_of_delta));
/*	for(i=0;i<size_of_delta;++i){
		printf("delta[%d]: %d \n", i, delta[i]);
	}
*/

	// Set d = H'(delta) and Concat c
	hash = calloc(LWE_N/8, sizeof(unsigned char));
	shake256(hash, LWE_N/8, delta, size_of_delta); // H'(delta) is LWE_N bit = LWE_N/8 bytes.
	for(i=0;i<LWE_N/8;++i){  // Concat C
		c[LWE_N+LWE_N+i]=hash[i];
	}
/*
	for(i=0;i<MESSAGE_LEN + LWE_N/8;++i){
		printf("hash[%d]: %d \n", i, hash[i]);
	}
*/


	// Set r = H(delta)
	unsigned char r[LWE_N];
	memset(r, 0, LWE_N);
	uint16_t r_idx[HR];
	unsigned int r_random_idx; 
	int hw, count = 0;

	hash = calloc(HR*4, sizeof(unsigned char));
	shake256(hash, HR*4, delta, size_of_delta); //같은 delta로 가능??

	while (hw < HR) {
		r_random_idx = hash[count++]; 
		r_random_idx <<= 8;	
		r_random_idx ^= hash[count++];
		r_random_idx = r_random_idx & (LWE_N - 1); // (seed1) || (seed2) 
		if (r[r_random_idx] == 0) {
			r[r_random_idx] = (hash[count++] & 0x02) - 1;
			hw++;
		}
		if (count >= (HR*4 - 3)) { 
			printf("New seed hash!\n");
			shake256(hash, LWE_N/8, hash, LWE_N/8);
			count = 0;
		}
	}
	//free(hash);
/*
	int one=0;
	int minus=0;
	int zero=0;
	for(i = 0; i<LWE_N; ++i){
		printf("r[%d]: %d \n", i, r[i]);
		if (r[i] == 0x01){
			one++;
		}
		else if (r[i] == 0xff){
			minus++;
		}
		else { zero++;}
	}
	printf("one: %d \t minus: %d \t zero: %d \n", one, minus, zero);
*/

	// Generate r_idx
	int neg_start = 0, back_position = HR;

	for (i = 0; i < LWE_N; ++i) {
		if (r[i] == 0x01){r_idx[neg_start++] = i;}
		else if (r[i] == 0xff){r_idx[--back_position] = i;}
	}
/*
	for(i = 0; i<HR; ++i){
		printf("r_idx[%d]: %d \n", i, r_idx[i]);
	}
*/


	// INSERT BCH CODE
	unsigned char delta_pad[size_of_delta + 1]; // delta_pad = delta||0x00, because BCH[511, 264, 59]
	delta_pad[size_of_delta]=0;
	for(i = 0; i < size_of_delta; ++i){
		delta_pad[i]=delta[i];
		//printf("delta_pad[%d]: %d \n", i, delta_pad[i]);
	}
	//printf("delta_pad[32]: %d \n", delta_pad[32]);

	unsigned char delta_hat[LWE_N / 8]; // delta_hat is result of encBCH(delta_pad)
	ecc_enc(delta_pad, delta_hat);
/*
	for(i = 0; i < LWE_N / 8; ++i){
		printf("delta_hat[%d]: %d \n", i, delta_hat[i]);
	}
*/
/*
	delta_hat[1]=0;
	delta_hat[2]=0;
	delta_hat[3]=0;

	for(i = 0; i < LWE_N / 8; ++i){
		printf("delta_hat[%d]: %d \n", i, delta_hat[i]);
	}
	unsigned char delta_pad_de[size_of_delta + 1];
	ecc_dec(delta_pad_de, delta_hat);

	for(i = 0; i < size_of_delta + 1; ++i){
		printf("delta_pad_de[%d]: %d \n", i, delta_pad_de[i]);
	}
*/


	// Parse pk_b & Make pk_a
	unsigned char pk_a[LWE_N];
	unsigned char pk_b[LWE_N];
	unsigned char seed_a[SEED_LEN];

	for(i=0; i<SEED_LEN; i++){
		seed_a[i] = pk[i];
	}
	shake256(pk_a, LWE_N, seed_a, SEED_LEN);

	for (i = 0; i < LWE_N; i++) {
		pk_b[i] = pk[SEED_LEN + i];
//		printf("pk_a[%d]: %d \t pk_b[%d]: %d \n", i, pk_a[i], i, pk_b[i]);
	}

	memset(c1h_b, 0, LWE_N*2);
	for (i = 0; i < LWE_N / 8; ++i) {
		for (j = 0; j < 8; ++j) {
			c1h_b[8 * i + j] = (delta_hat[i] >> j) << _8_LOG_T;// Shift each bit of delta_hat to MSB and save in c2h_b
		}
	}
/*
	for (i = 0; i < LWE_N*2; i++) {
		printf("c2h_b[%d]: %d \n", i, c2h_b[i]);
	}
*/

	// Compute a * r and b * r, and then add to c2h_a and c2h_b, respectively.


	//printf("neg_start: %d \n", neg_start);

	memset(c1h_a, 0, LWE_N*2);
	for(i = 0; i < HR; ++i){
		uint16_t branch = (2 * ((i - neg_start) >> sft & 0x1) - 1);
		uint16_t deg = r_idx[i];
		for(j = 0; j < LWE_N; ++j){
			c1h_a[deg+j] += branch * pk_a[j];
			c1h_b[deg+j] += branch * pk_b[j];
		}
	}
	for(j = 0; j < LWE_N; ++j){
		c1h_a[j] -= c1h_a[LWE_N+j];
		c1h_b[j] -= c1h_b[LWE_N+j];
	}
/*
	for(j = 0; j < LWE_N; ++j){
		printf("c2h_a[%d]: %d \t c2h_b[%d]: %d \n", j, c2h_a[j], j, c2h_b[j]);
	}
*/
	// Send c2h_a and c2h_b from mod q to mod p and mod k

	for (i=0; i< LWE_N; ++i) {
		c[i] = ((c1h_a[i] + 0x02) & 0xfc);
		c[LWE_N + i] = ((c1h_b[i] + 0x04) & 0xf8);
	}

	// G(c1, d, delta)
	hash_t=calloc((LWE_N+LWE_N+LWE_N/8+size_of_delta), sizeof(unsigned char));
	memcpy(hash_t, c, (LWE_N+LWE_N+LWE_N/8));
	memcpy(hash_t+(LWE_N+LWE_N+LWE_N/8), delta, size_of_delta);
	shake256(shared_k, MESSAGE_LEN, hash_t, LWE_N+LWE_N+LWE_N/8+size_of_delta);
/*
	for (i=0;i<LWE_N;++i){
		printf("c1h_a[%d]: %d \t c2h_b[%d]: %d \n", i, c[i], i, c[LWE_N+i]);
	}
	for (i=0;i<LWE_N/8;++i){
		printf("d[%d]: %d \n", i, c[LWE_N+LWE_N+i]);
	}
	for (i=0;i<size_of_delta;++i){
		printf("delta[%d]: %d \n", i, delta[i]);
	}
*/


/*
	for (i=0;i<MESSAGE_LEN+LWE_N+LWE_N+LWE_N/8;++i){
		printf("c[%d]: %d \n", i, c[i]);
	}
*/
/*
	for (i=MESSAGE_LEN;i<MESSAGE_LEN+LWE_N+LWE_N;++i){
		printf("c[%d]: %d \n", i, c[i]);
	}
*/
	free(hash);
	free(hash_t);

	finish = clock();
	elapsed1 += (finish - start);
	return 0;
}




int Dec(unsigned char *shared_k, unsigned char *c, const unsigned char *sk, const unsigned char *pk){
	start = clock();
	int res = 0;
	int i, j;
	//unsigned int tmp;
	unsigned char r[LWE_N] = { 0, };
	unsigned int r_idx[HR];

	//uint64_t delta[LWE_N / 64] = { 0, };
	unsigned char delta_hat[LWE_N/8] = { 0, };
	unsigned char delta_pad[size_of_delta + 1] = {0,};
	unsigned char delta[size_of_delta] = { 0, };
	unsigned char *hash = NULL;
	unsigned char *hash_t = NULL;

	// Initialize c1h_b as q/2 * delta_hat
	unsigned char c1h_a[LWE_N*2] = { 0, };
	unsigned char c1h_b[LWE_N*2] = { 0, };
	unsigned char decomp_delta[LWE_N*2]={0,};

	// Initialize delta as c1h_b
	memset(decomp_delta, 0, LWE_N*2);
	memset(c1h_a, 0, LWE_N*2);
	for (i=0; i<LWE_N; ++i) {
		decomp_delta[i] = c[LWE_N + i];
		c1h_a[i] = c[i];
	}
/*
	for(i=0;i<LWE_N*2;++i){
		printf("c2h_a[%d]: %d \n", i, c2h_a[i]);
	}
*/

/*
	for(i=0;i<LWE_N*2;++i){
		printf("decomp_delta(c2h_b)[%d]: %d \t c2h_a[%d]: %d \n", i, decomp_delta[i], i, c2h_a[i]);
	}
*/
	//DEC algorithm step2 NONE....


  // gen s_idx
	uint16_t sk_s[HS];
	int neg_start = 0, back_position = HS;	
	for (i = 0; i < LWE_N; ++i) {
		if (sk[i] == 0x01){ 
			sk_s[neg_start++] = i;
		}
		if (sk[i] == 0xff){
			sk_s[--back_position] = i;
		}
	}

	//gen u (parsing sk)
	unsigned char u[LWE_N/8]={0,};
	for (i=LWE_N;i<LWE_N+LWE_N/8;++i){
		u[i-LWE_N]=sk[i];
	}

/*
	for(i = 0; i<HS; ++i){
		printf("sk_s[%d]: %d \n", i, sk_s[i]);
	}
	for(i = 0; i<LWE_N; ++i){
		printf("sk[%d]: %d \n", i, sk[i]);
	}
*/

	// Compute delta (delta(= c1b) + c1(=c1a) * s)
	for(i = 0; i < HS; ++i){
		uint16_t branch = (2 * ((i - neg_start) >> sft & 0x1) - 1);
		uint16_t deg = sk_s[i];
		for(int j = 0; j < LWE_N; ++j){
			decomp_delta[deg+j] += branch * c1h_a[j];
	    }
	}
	for(j = 0; j < LWE_N; ++j){
		decomp_delta[j] -= decomp_delta[LWE_N+j];
	}
/*
	for(j = 0; j < LWE_N; ++j){
		printf("decomp_delta[%d]: %d \n", j, decomp_delta[j]);
	}
*/
	// Compute delta = 2/p * delta
	for (i = 0; i < LWE_N; ++i) {
		decomp_delta[i] += 0x40;
		decomp_delta[i] >>= _8_LOG_T;
	}
	// Set delta_hat
	for (i = 0; i < LWE_N/8; ++i) {
		for (j = 0; j < 8; ++j) {
			uint8_t a = (decomp_delta[8 * i + j]) << j;
			delta_hat[i] ^= a;  
		}
	}
/*
	for(i=0;i<LWE_N/8;++i){
		printf("delta_hat[%d]: %d \n", i, delta_hat[i]);
	}
*/

	// INSERT BCH DECODE ///////////////////////////////////////////////
	ecc_dec(delta_pad, delta_hat);
/*
	for (i=0;i<size_of_delta+1;++i){
		printf("delta_pad[%d]: %d \n", i, delta_pad[i]);
	}
*/

	for (i = 0; i < size_of_delta; ++i) {
		delta[i] = delta_pad[i];
	}

	// Set d = H'(delta)
	hash = calloc(LWE_N/8, sizeof(unsigned char));
	unsigned char d[LWE_N/8]={0,};
	shake256(hash, LWE_N/8, delta, size_of_delta); // H'(delta) is LWE_N bit = LWE_N/8 bytes.
	for(i=0;i<LWE_N/8;++i){  
		d[i]=hash[i];
	}

	// Set r = H(delta)
	memset(r, 0, LWE_N);
	unsigned int r_random_idx; 
	int hw, count = 0;

	hash = calloc(HR*4, sizeof(unsigned char));
	shake256(hash, HR*4, delta, size_of_delta); //같은 delta로 가능??

	while (hw < HR) {
		r_random_idx = hash[count++]; 
		r_random_idx <<= 8;	
		r_random_idx ^= hash[count++];
		r_random_idx = r_random_idx & (LWE_N - 1); // (seed1) || (seed2) 
		if (r[r_random_idx] == 0) {
			r[r_random_idx] = (hash[count++] & 0x02) - 1;
			hw++;
		}
		if (count >= (HR*4 - 3)) { 
			printf("DEC_New seed hash!\n");
			shake256(hash, LWE_N/8, hash, LWE_N/8);
			count = 0;
		}
	}
	
/*
	int one=0;
	int minus=0;
	int zero=0;
	for(i = 0; i<LWE_N; ++i){
		printf("r[%d]: %d \n", i, r[i]);
		if (r[i] == 0x01){
			one++;
		}
		else if (r[i] == 0xff){
			minus++;
		}
		else { zero++;}
	}
	printf("one: %d \t minus: %d \t zero: %d \n", one, minus, zero);
*/


	// Generate r_idx
	neg_start = 0;
	back_position = HR;

	for (i = 0; i < LWE_N; ++i) {
		if (r[i] == 0x01){r_idx[neg_start++] = i;}
		else if (r[i] == 0xff){r_idx[--back_position] = i;}
	}

/*
	for(i = 0; i<HR; ++i){
		printf("r_idx[%d]: %d \n", i, r_idx[i]);
	}
*/

	// INSERT BCH CODE
	memset(delta_pad, 0, size_of_delta+1);
	for(i = 0; i < size_of_delta; ++i){
		delta_pad[i]=delta[i];
		//printf("delta_pad[%d]: %d \n", i, delta_pad[i]);
	}

	memset(delta_hat, 0, LWE_N/8);
	ecc_enc(delta_pad, delta_hat);


	// Parse pk_b & Make pk_a
	unsigned char pk_a[LWE_N];
	unsigned char pk_b[LWE_N];
	unsigned char seed_a[SEED_LEN];

	for(i=0; i<SEED_LEN; i++){
		seed_a[i] = pk[i];
	}
	shake256(pk_a, LWE_N, seed_a, SEED_LEN);

	for (i = 0; i < LWE_N; i++) {
		pk_b[i] = pk[SEED_LEN + i];
//		printf("pk_a[%d]: %d \t pk_b[%d]: %d \n", i, pk_a[i], i, pk_b[i]);
	}
	memset(c1h_b, 0, LWE_N*2);
	for (i = 0; i < LWE_N / 8; ++i) {
		for (j = 0; j < 8; ++j) {
			c1h_b[8 * i + j] = (delta_hat[i] >> j) << _8_LOG_T;// Shift each bit of delta_hat to MSB and save in c2h_b
		}
	}
/*
	for (i = 0; i < LWE_N*2; i++) {
		printf("c2h_b[%d]: %d \n", i, c2h_b[i]);
	}
*/

	// Compute a * r and b * r, and then add to c2h_a and c2h_b, respectively.


	//printf("neg_start: %d \n", neg_start);

	memset(c1h_a, 0, LWE_N*2);
	for(i = 0; i < HR; ++i){
		uint16_t branch = (2 * ((i - neg_start) >> sft & 0x1) - 1);
		uint16_t deg = r_idx[i];
		for(j = 0; j < LWE_N; ++j){
			c1h_a[deg+j] += branch * pk_a[j];
			c1h_b[deg+j] += branch * pk_b[j];
		}
	}
	for(j = 0; j < LWE_N; ++j){
		c1h_a[j] -= c1h_a[LWE_N+j];
		c1h_b[j] -= c1h_b[LWE_N+j];
	}
/*
	for(j = 0; j < LWE_N; ++j){
		printf("c2h_a[%d]: %d \t c2h_b[%d]: %d \n", j, c2h_a[j], j, c2h_b[j]);
	}
*/
	// Send c2h_a and c2h_b from mod q to mod p and mod k


	for (i = 0; i < LWE_N; ++i) {
		c1h_a[i] = (c1h_a[i] + 0x02) & 0xfc;
		c1h_b[i] = (c1h_b[i] + 0x04) & 0xf8;
	}

	for (i = 0; i < LWE_N; ++i) {
		if ((c[i] != ((c1h_a[i] + 0x02) & 0xfc)) || (c[LWE_N + i] != ((c1h_b[i] + 0x04) & 0xf8))){
			printf("c1h_a or c1h_b error!!!!");		

			// G(c1, d, u)
			hash_t=calloc((LWE_N+LWE_N+LWE_N/8+LWE_N/8), sizeof(unsigned char));
			memcpy(hash_t, c1h_a, LWE_N);
			memcpy(hash_t+LWE_N, c1h_b, LWE_N);
			memcpy(hash_t+LWE_N+LWE_N, d, LWE_N/8);
			memcpy(hash_t+LWE_N+LWE_N+LWE_N/8, u, LWE_N/8);

			shake256(shared_k, MESSAGE_LEN, hash_t, LWE_N+LWE_N+LWE_N/8+LWE_N/8);

			free(hash_t);
	
			free(hash);
			return res = 2;
		}
	}

	//G(c1, d, delta)
	hash_t=calloc((LWE_N+LWE_N+LWE_N/8+size_of_delta), sizeof(unsigned char));
	memcpy(hash_t, c1h_a, LWE_N);
	memcpy(hash_t+LWE_N, c1h_b, LWE_N);
	memcpy(hash_t+LWE_N+LWE_N, d, LWE_N/8);
	memcpy(hash_t+LWE_N+LWE_N+LWE_N/8, delta, size_of_delta);
	shake256(shared_k, MESSAGE_LEN, hash_t, LWE_N+LWE_N+LWE_N/8+size_of_delta);
/*
	for (i=0;i<LWE_N;++i){
		printf("c1h_a[%d]: %d \t c2h_b[%d]: %d \n", i, c1h_a[i], i, c1h_b[i]);
	}
	for (i=0;i<LWE_N/8;++i){
		printf("d[%d]: %d \n", i, d[i]);
	}
	for (i=0;i<size_of_delta;++i){
		printf("delta[%d]: %d \n", i, delta[i]);
	}
*/


	free(hash_t);


	free(hash);

	finish = clock();
	elapsed2 += (finish - start);
	return res;
}
