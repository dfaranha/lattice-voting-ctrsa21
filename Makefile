CPP = g++
CFLAGS = -O3 -march=native -mtune=native -Wall -ggdb -pthread
INCLUDES = bench.h cpucycles.h
BENCH = bench.c cpucycles.c
TEST = test.c
GAUSSIAN = gaussian.o fastrandombytes.c randombytes.c
LIBS = -lflint -lgmp

.PHONY: all test bench clean

all: commit encrypt vericrypt shuffle

# Run only the test phase of each binary. The benchmarks dominate the runtime,
# so keeping them out of a test run is what makes this usable in CI.
test: all
	./commit test
	./encrypt test
	./vericrypt test
	./shuffle test

bench: all
	./commit bench
	./encrypt bench
	./vericrypt bench
	./shuffle bench

commit: commit.c ${TEST} ${BENCH} ${INCLUDES} gaussian_ct.cpp
	${CPP} ${CFLAGS} -DSIGMA_PARAM=SIGMA_C -c gaussian_ct.cpp -o gaussian.o
	${CPP} ${CFLAGS} -DMAIN commit.c ${GAUSSIAN} ${TEST} ${BENCH} -o commit ${LIBS}

encrypt: encrypt.c ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -DMAIN encrypt.c ${TEST} ${BENCH} -o encrypt ${LIBS}

vericrypt: vericrypt.c encrypt.c ${TEST} ${BENCH} ${INCLUDES} gaussian_ct.cpp
	${CPP} ${CFLAGS} -DSIGMA_PARAM=SIGMA_E -c gaussian_ct.cpp -o gaussian.o
	${CPP} ${CFLAGS} -c encrypt.c -o encrypt.o
	${CPP} ${CFLAGS} -DMAIN vericrypt.c encrypt.o sha224-256.c ${GAUSSIAN} ${TEST} ${BENCH} -o vericrypt ${LIBS}

shuffle: shuffle.c commit.c ${TEST} ${BENCH}
	${CPP} ${CFLAGS} -DSIGMA_PARAM=SIGMA_C -c gaussian_ct.cpp -o gaussian.o
	${CPP} ${CFLAGS} commit.c shuffle.c sha224-256.c ${GAUSSIAN} ${TEST} ${BENCH} -o shuffle ${LIBS}

clean:
	rm -f *.o commit encrypt vericrypt shuffle
