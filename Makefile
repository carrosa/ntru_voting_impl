CPP = g++ -g3
CFLAGS = -O3 -march=native -mtune=native -Wall -ggdb -I NFLlib/include/ -I NFLlib/include/nfl -I NFLlib/include/nfl/prng -I include -DNFL_OPTIMIZED=ON -DNTT_AVX2 -flto
INCLUDES = include/bench.h include/cpucycles.h
BLAKE3 = src/blake3/blake3.c src/blake3/blake3_dispatch.c src/blake3/blake3_portable.c src/blake3/blake3_sse2_x86-64_unix.S src/blake3/blake3_sse41_x86-64_unix.S src/blake3/blake3_avx2_x86-64_unix.S src/blake3/blake3_avx512_x86-64_unix.S
BENCH = src/bench.c src/cpucycles.c
TEST = src/test.c
LIBS = deps/libnfllib_static.a -lgmp -lmpfr -L deps/ -lflint -lquadmath

all: bdlop ntru_bdlop ntru ntru_shuffle ntru_pismall

bdlop: src/bdlop.cpp src/bgv.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -c src/bgv.cpp -o bgv.o
	${CPP} ${CFLAGS} -DMAIN src/bdlop.cpp bgv.o ${TEST} ${BENCH} -o bdlop ${LIBS}

ntru_bdlop: src/ntru_bdlop.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -DMAIN src/ntru_bdlop.cpp ${TEST} ${BENCH} -o ntru_bdlop ${LIBS}

ntru: src/ntru.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -c src/sample_z_small.c -o sample_z_small.o
	${CPP} ${CFLAGS} -DMAIN src/ntru.cpp sample_z_small.o ${TEST} ${BENCH} ${BLAKE3} -o ntru ${LIBS}

ntru_pismall: src/bdlop.cpp src/ntru_pismall.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -DSIZE=3 -c src/bdlop.cpp -o ntru_bdlop.o
	${CPP} ${CFLAGS} -DSIZE=3 -DMAIN src/ntru_pismall.cpp ntru_bdlop.o ${TEST} ${BENCH} ${BLAKE3} -o ntru_pismall ${LIBS}

ntru_shuffle: src/ntru_shuffle.cpp src/ntru_bdlop.cpp ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -c src/sample_z_small.c -o sample_z_small.o
	${CPP} ${CFLAGS} -c src/ntru_bdlop.cpp -o ntru_bdlop.o
	${CPP} ${CFLAGS} -DMAIN src/ntru_shuffle.cpp sample_z_small.o ntru_bdlop.o ${TEST} ${BENCH} ${BLAKE3} -o ntru_shuffle ${LIBS}

clean:
	rm -f *.o bdlop ntru ntru_shuffle ntru_pismall ntru_bdlop
