CPP = g++
# Build-time constants, e.g. make CONFIG=-DMSGS=8 voting for an election of 8
# voters, or CONFIG=-DTESTS=5 for fewer test iterations.
CONFIG =
CFLAGS = -O3 -march=native -mtune=native -Wall -ggdb -pthread $(CONFIG)
INCLUDES = bench.h cpucycles.h param.h
BENCH = bench.c cpucycles.c
TEST = test.c
RAND = fastrandombytes.c randombytes.c
LIBS = -lflint -lgmp

# LaZer, for the voters' validity proofs (validity.c), is the submodule
# third_party/lazer: branch skip-zeros of dfaranha/lazer, LaZer with fixes and
# speed-ups. Its own Makefile builds it on first use, HEXL included; fetch it
# with git submodule update --init third_party/lazer. validity.c is C, compiled
# with CC against LaZer's headers; the rest, built as C++, calls it through
# validity.h.
CC = gcc
LAZER = third_party/lazer
LAZER_LIB = $(LAZER)/liblazer.a
HEXL_LIB = $(LAZER)/third_party/hexl-development/build/hexl/lib64/libhexl.a
LAZER_LIBS = $(LAZER_LIB) $(HEXL_LIB) -lmpfr -lstdc++

.PHONY: all test bench clean FORCE

all: commit encrypt vericrypt shuffle voting

# make compares timestamps, so it cannot see that CONFIG changed since the
# binaries were built: record CONFIG in a stamp everything depends on,
# rewritten only when it changes.
STAMP = .config
$(STAMP): FORCE
	@printf '%s' '$(CONFIG)' | cmp -s - $@ 2>/dev/null || printf '%s' '$(CONFIG)' > $@

# Run only the test phase of each binary. The benchmarks dominate the runtime,
# so keeping them out of a test run is what makes this usable in CI. TESTS,
# BENCH and MSGS can be overridden at build time for a faster configuration,
# for example CONFIG=-DTESTS=5.
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
gaussian_c.o: gaussian_ct.cpp fastrandombytes.h ${INCLUDES} $(STAMP)
	${CPP} ${CFLAGS} -DSIGMA_PARAM=SIGMA_C -c gaussian_ct.cpp -o $@

gaussian_e.o: gaussian_ct.cpp fastrandombytes.h ${INCLUDES} $(STAMP)
	${CPP} ${CFLAGS} -DSIGMA_PARAM=SIGMA_E -c gaussian_ct.cpp -o $@

# The proof of shuffle also masks the committed permutation elements, with a
# much narrower Gaussian. Renaming the entry point lets both widths be linked
# into the same binary.
gaussian_s.o: gaussian_ct.cpp fastrandombytes.h ${INCLUDES} $(STAMP)
	${CPP} ${CFLAGS} -DSIGMA_PARAM=SIGMA_S \
		-Ddiscrete_gaussian=discrete_gaussian_small -c gaussian_ct.cpp -o $@

# The end-to-end test links the verifiable encryption, at SIGMA_E, with the
# commitments and the proof of shuffle, at SIGMA_C: both samplers are called
# discrete_gaussian, so the verifiable encryption gets one of its own name.
gaussian_ev.o: gaussian_ct.cpp fastrandombytes.h ${INCLUDES} $(STAMP)
	${CPP} ${CFLAGS} -DSIGMA_PARAM=SIGMA_E \
		-Ddiscrete_gaussian=discrete_gaussian_e -c gaussian_ct.cpp -o $@

vericrypt_v.o: vericrypt.c vericrypt.h encrypt.h commit.h ${INCLUDES} $(STAMP)
	${CPP} ${CFLAGS} -Ddiscrete_gaussian=discrete_gaussian_e -c vericrypt.c -o $@

encrypt.o: encrypt.c encrypt.h ${INCLUDES} $(STAMP)
	${CPP} ${CFLAGS} -c encrypt.c -o $@

commit: commit.c commit.h sha224-256.c ${RAND} fastrandombytes.h ${TEST} ${BENCH} ${INCLUDES} gaussian_c.o gaussian_s.o $(STAMP)
	${CPP} ${CFLAGS} -DMAIN commit.c sha224-256.c gaussian_c.o gaussian_s.o ${RAND} ${TEST} ${BENCH} -o $@ ${LIBS}

encrypt: encrypt.c encrypt.h ${TEST} ${BENCH} ${INCLUDES} $(STAMP)
	${CPP} ${CFLAGS} -DMAIN encrypt.c ${TEST} ${BENCH} -o $@ ${LIBS}

vericrypt: vericrypt.c encrypt.o sha224-256.c ${RAND} fastrandombytes.h ${TEST} ${BENCH} ${INCLUDES} gaussian_e.o $(STAMP)
	${CPP} ${CFLAGS} -DMAIN vericrypt.c encrypt.o sha224-256.c gaussian_e.o ${RAND} ${TEST} ${BENCH} -o $@ ${LIBS}

shuffle: shuffle.c commit.c commit.h serial.c serial.h sha224-256.c ${RAND} fastrandombytes.h ${TEST} ${BENCH} ${INCLUDES} gaussian_c.o gaussian_s.o $(STAMP)
	${CPP} ${CFLAGS} commit.c serial.c shuffle.c sha224-256.c gaussian_c.o gaussian_s.o ${RAND} ${TEST} ${BENCH} -o $@ ${LIBS}

# The voting protocol end to end: voters, board, authority and verifier.
$(LAZER_LIB):
	$(MAKE) -C $(LAZER) liblazer.a
$(HEXL_LIB): $(LAZER_LIB)

# The validity proof, in C against LaZer (its parameters: lazer/validity.h,
# generated from lazer/validity.py by LaZer's lin-codegen.sage).
validity.o: validity.c validity.h lazer/validity.h $(LAZER_LIB)
	$(CC) -O3 -march=native -mtune=native -Wall -Wextra -ggdb -I $(LAZER) -c validity.c -o $@

voting: voting.c shuffle.c commit.c commit.h serial.c serial.h vericrypt.h encrypt.h validity.h vericrypt_v.o encrypt.o gaussian_ev.o validity.o sha224-256.c ${RAND} fastrandombytes.h ${TEST} ${BENCH} ${INCLUDES} gaussian_c.o gaussian_s.o $(STAMP)
	${CPP} ${CFLAGS} voting.c commit.c serial.c sha224-256.c vericrypt_v.o encrypt.o gaussian_c.o gaussian_s.o gaussian_ev.o validity.o ${RAND} ${TEST} ${BENCH} -o $@ ${LIBS} $(LAZER_LIBS)

# Fuzzing the proof decoder, with clang's libFuzzer under AddressSanitizer and
# UndefinedBehaviorSanitizer; fuzz_serial.c says how to run it.
FUZZ_CPP = clang++
FUZZ_FLAGS = -O1 -g -march=native -pthread -fsanitize=address,undefined \
	-fno-sanitize-recover=all

# The samplers are GCC-only and the decoder does not call them, so their
# objects come from the ordinary build.
fuzz_serial: fuzz_serial.c serial.c serial.h commit.c commit.h sha224-256.c ${RAND} fastrandombytes.h ${INCLUDES} gaussian_c.o gaussian_s.o $(STAMP)
	${FUZZ_CPP} ${FUZZ_FLAGS} -fsanitize=fuzzer -x c++ fuzz_serial.c serial.c commit.c sha224-256.c ${RAND} -x none gaussian_c.o gaussian_s.o -o $@ ${LIBS}

clean:
	rm -f *.o commit encrypt vericrypt shuffle voting fuzz_serial $(STAMP)
