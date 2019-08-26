CC = gcc

CFLAGS=-O3 -fomit-frame-pointer -msse2avx -mavx2 -march=native -std=c99

all :
	$(CC) $(CFLAGS) -c RLizarMon_KEM_Comfort.c main.c randombytes.c fips202.c bch.c ecc.c
	$(CC) $(CFLAGS) -o RLizarMon_KEM_Comfort RLizarMon_KEM_Comfort.o main.o randombytes.o fips202.o bch.o ecc.o
	
run : all
	./RLizarMon_KEM_Comfort

clean :
	rm -f *.o
	rm -f RLizarMon_KEM_Comfort

new :
	make clean
	make all
	./RLizarMon_KEM_Comfort
