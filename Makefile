CPP = g++
CFLAGS = -O3 -march=native -mtune=native -Wall -ggdb -pthread
INCLUDES = bench.h cpucycles.h param.h
BENCH = bench.c cpucycles.c
TEST = test.c
RAND = fastrandombytes.c randombytes.c
LIBS = -lflint -lgmp

.PHONY: all test bench clean

all: commit encrypt vericrypt shuffle voting

# Run only the test phase of each binary. The benchmarks dominate the runtime,
# so keeping them out of a test run is what makes this usable in CI. TESTS,
# BENCH and MSGS can be overridden at build time for a faster configuration,
# for example CFLAGS="... -DTESTS=5".
test: all
	./commit test
	./encrypt test
	./vericrypt test
	./shuffle test
	./voting

bench: all
	./commit bench
	./encrypt bench
	./vericrypt bench
	./shuffle bench

# The discrete Gaussian sampler bakes its standard deviation in at compile
# time, so each width needs its own object file. Sharing a single gaussian.o
# between targets made a parallel build race on that file, and could link a
# sampler of the wrong width into a binary.
gaussian_c.o: gaussian_ct.cpp fastrandombytes.h ${INCLUDES}
	${CPP} ${CFLAGS} -DSIGMA_PARAM=SIGMA_C -c gaussian_ct.cpp -o $@

gaussian_e.o: gaussian_ct.cpp fastrandombytes.h ${INCLUDES}
	${CPP} ${CFLAGS} -DSIGMA_PARAM=SIGMA_E -c gaussian_ct.cpp -o $@

# The proof of shuffle also masks the committed permutation elements, with a
# much narrower Gaussian. Renaming the entry point lets both widths be linked
# into the same binary.
gaussian_s.o: gaussian_ct.cpp fastrandombytes.h ${INCLUDES}
	${CPP} ${CFLAGS} -DSIGMA_PARAM=SIGMA_S \
		-Ddiscrete_gaussian=discrete_gaussian_small -c gaussian_ct.cpp -o $@

# The end-to-end test links the verifiable encryption, at SIGMA_E, with the
# commitments and the proof of shuffle, at SIGMA_C: both samplers are called
# discrete_gaussian, so the verifiable encryption gets one of its own name.
gaussian_ev.o: gaussian_ct.cpp fastrandombytes.h ${INCLUDES}
	${CPP} ${CFLAGS} -DSIGMA_PARAM=SIGMA_E \
		-Ddiscrete_gaussian=discrete_gaussian_e -c gaussian_ct.cpp -o $@

vericrypt_v.o: vericrypt.c vericrypt.h encrypt.h commit.h ${INCLUDES}
	${CPP} ${CFLAGS} -Ddiscrete_gaussian=discrete_gaussian_e -c vericrypt.c -o $@

encrypt.o: encrypt.c encrypt.h ${INCLUDES}
	${CPP} ${CFLAGS} -c encrypt.c -o $@

commit: commit.c commit.h sha224-256.c ${RAND} fastrandombytes.h ${TEST} ${BENCH} ${INCLUDES} gaussian_c.o gaussian_s.o
	${CPP} ${CFLAGS} -DMAIN commit.c sha224-256.c gaussian_c.o gaussian_s.o ${RAND} ${TEST} ${BENCH} -o $@ ${LIBS}

encrypt: encrypt.c encrypt.h ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -DMAIN encrypt.c ${TEST} ${BENCH} -o $@ ${LIBS}

vericrypt: vericrypt.c encrypt.o sha224-256.c ${RAND} fastrandombytes.h ${TEST} ${BENCH} ${INCLUDES} gaussian_e.o
	${CPP} ${CFLAGS} -DMAIN vericrypt.c encrypt.o sha224-256.c gaussian_e.o ${RAND} ${TEST} ${BENCH} -o $@ ${LIBS}

shuffle: shuffle.c commit.c commit.h serial.c serial.h sha224-256.c ${RAND} fastrandombytes.h ${TEST} ${BENCH} ${INCLUDES} gaussian_c.o gaussian_s.o
	${CPP} ${CFLAGS} commit.c serial.c shuffle.c sha224-256.c gaussian_c.o gaussian_s.o ${RAND} ${TEST} ${BENCH} -o $@ ${LIBS}

# The voting protocol end to end: voters, board, authority and verifier.
voting: voting.c shuffle.c commit.c commit.h serial.c serial.h vericrypt.h encrypt.h vericrypt_v.o encrypt.o gaussian_ev.o sha224-256.c ${RAND} fastrandombytes.h ${TEST} ${BENCH} ${INCLUDES} gaussian_c.o gaussian_s.o
	${CPP} ${CFLAGS} voting.c commit.c serial.c sha224-256.c vericrypt_v.o encrypt.o gaussian_c.o gaussian_s.o gaussian_ev.o ${RAND} ${TEST} ${BENCH} -o $@ ${LIBS}

# Fuzzing the proof decoder, with clang's libFuzzer under AddressSanitizer and
# UndefinedBehaviorSanitizer; fuzz_serial.c says how to run it.
FUZZ_CPP = clang++
FUZZ_FLAGS = -O1 -g -march=native -pthread -fsanitize=address,undefined \
	-fno-sanitize-recover=all

# The samplers are GCC-only and the decoder does not call them, so their
# objects come from the ordinary build.
fuzz_serial: fuzz_serial.c serial.c serial.h commit.c commit.h sha224-256.c ${RAND} fastrandombytes.h ${INCLUDES} gaussian_c.o gaussian_s.o
	${FUZZ_CPP} ${FUZZ_FLAGS} -fsanitize=fuzzer -x c++ fuzz_serial.c serial.c commit.c sha224-256.c ${RAND} -x none gaussian_c.o gaussian_s.o -o $@ ${LIBS}

clean:
	rm -f *.o commit encrypt vericrypt shuffle voting fuzz_serial
