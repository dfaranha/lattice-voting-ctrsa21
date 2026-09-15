CPP = g++
CFLAGS = -O3 -march=native -mtune=native -Wall -ggdb -pthread
INCLUDES = bench.h cpucycles.h param.h
BENCH = bench.c cpucycles.c
TEST = test.c
RAND = fastrandombytes.c randombytes.c
LIBS = -lflint -lgmp

.PHONY: all test bench clean

all: commit encrypt vericrypt shuffle

# Run only the test phase of each binary. The benchmarks dominate the runtime,
# so keeping them out of a test run is what makes this usable in CI. TESTS,
# BENCH and MSGS can be overridden at build time for a faster configuration,
# for example CFLAGS="... -DTESTS=5".
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

# The discrete Gaussian sampler bakes its standard deviation in at compile
# time, so each width needs its own object file. Sharing a single gaussian.o
# between targets made a parallel build race on that file, and could link a
# sampler of the wrong width into a binary.
gaussian_c.o: gaussian_ct.cpp ${INCLUDES}
	${CPP} ${CFLAGS} -DSIGMA_PARAM=SIGMA_C -c gaussian_ct.cpp -o $@

gaussian_e.o: gaussian_ct.cpp ${INCLUDES}
	${CPP} ${CFLAGS} -DSIGMA_PARAM=SIGMA_E -c gaussian_ct.cpp -o $@

encrypt.o: encrypt.c encrypt.h ${INCLUDES}
	${CPP} ${CFLAGS} -c encrypt.c -o $@

commit: commit.c commit.h ${TEST} ${BENCH} ${INCLUDES} gaussian_c.o
	${CPP} ${CFLAGS} -DMAIN commit.c gaussian_c.o ${RAND} ${TEST} ${BENCH} -o $@ ${LIBS}

encrypt: encrypt.c encrypt.h ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -DMAIN encrypt.c ${TEST} ${BENCH} -o $@ ${LIBS}

vericrypt: vericrypt.c encrypt.o ${TEST} ${BENCH} ${INCLUDES} gaussian_e.o
	${CPP} ${CFLAGS} -DMAIN vericrypt.c encrypt.o sha224-256.c gaussian_e.o ${RAND} ${TEST} ${BENCH} -o $@ ${LIBS}

shuffle: shuffle.c commit.c commit.h serial.c serial.h ${TEST} ${BENCH} ${INCLUDES} gaussian_c.o
	${CPP} ${CFLAGS} commit.c serial.c shuffle.c sha224-256.c gaussian_c.o ${RAND} ${TEST} ${BENCH} -o $@ ${LIBS}

clean:
	rm -f *.o commit encrypt vericrypt shuffle
