CPP = g++
CFLAGS = -O3 -march=native -mtune=native -Wall -ggdb -pthread
INCLUDES = bench.h cpucycles.h param.h
BENCH = bench.c cpucycles.c
TEST = test.c
RAND = fastrandombytes.c randombytes.c
LIBS = -lflint -lgmp

.PHONY: all test bench clean

all: commit encrypt vericrypt shuffle lnp

# Run only the test phase of each binary. The benchmarks dominate the runtime,
# so keeping them out of a test run is what makes this usable in CI. TESTS,
# BENCH and MSGS can be overridden at build time for a faster configuration,
# for example CFLAGS="... -DTESTS=5".
test: all
	./commit test
	./encrypt test
	./vericrypt test
	./shuffle test
	./lnp test

bench: all
	./commit bench
	./encrypt bench
	./vericrypt bench
	./shuffle bench
	./lnp bench

# The discrete Gaussian sampler bakes its standard deviation in at compile
# time, so each width needs its own object file. Sharing a single gaussian.o
# between targets made a parallel build race on that file, and could link a
# sampler of the wrong width into a binary.
gaussian_c.o: gaussian_ct.cpp ${INCLUDES}
	${CPP} ${CFLAGS} -DSIGMA_PARAM=SIGMA_C -c gaussian_ct.cpp -o $@

gaussian_e.o: gaussian_ct.cpp ${INCLUDES}
	${CPP} ${CFLAGS} -DSIGMA_PARAM=SIGMA_E -c gaussian_ct.cpp -o $@

# The proof of shuffle also masks the committed permutation elements, with a
# much narrower Gaussian. Renaming the entry point lets both widths be linked
# into the same binary.
gaussian_s.o: gaussian_ct.cpp ${INCLUDES}
	${CPP} ${CFLAGS} -DSIGMA_PARAM=SIGMA_S \
		-Ddiscrete_gaussian=discrete_gaussian_small -c gaussian_ct.cpp -o $@

# The approximate range proof masks a 256-coordinate projection, which wants a
# width between the two above.
gaussian_p.o: gaussian_ct.cpp ${INCLUDES}
	${CPP} ${CFLAGS} -DSIGMA_PARAM=SIGMA_P \
		-Ddiscrete_gaussian=discrete_gaussian_proj -c gaussian_ct.cpp -o $@

# The batched proof shares one challenge across all messages, so its masks are
# wider; see SIGMA_B in param.h.
gaussian_b.o: gaussian_ct.cpp ${INCLUDES}
	${CPP} ${CFLAGS} -DSIGMA_PARAM=SIGMA_B \
		-Ddiscrete_gaussian=discrete_gaussian_batch -c gaussian_ct.cpp -o $@

encrypt.o: encrypt.c encrypt.h ${INCLUDES}
	${CPP} ${CFLAGS} -c encrypt.c -o $@

commit: commit.c commit.h ${TEST} ${BENCH} ${INCLUDES} gaussian_c.o gaussian_s.o
	${CPP} ${CFLAGS} -DMAIN commit.c gaussian_c.o gaussian_s.o ${RAND} ${TEST} ${BENCH} -o $@ ${LIBS}

encrypt: encrypt.c encrypt.h ${TEST} ${BENCH} ${INCLUDES}
	${CPP} ${CFLAGS} -DMAIN encrypt.c ${TEST} ${BENCH} -o $@ ${LIBS}

vericrypt: vericrypt.c encrypt.o ${TEST} ${BENCH} ${INCLUDES} gaussian_e.o
	${CPP} ${CFLAGS} -DMAIN vericrypt.c encrypt.o sha224-256.c gaussian_e.o ${RAND} ${TEST} ${BENCH} -o $@ ${LIBS}

shuffle: shuffle.c commit.c commit.h lnp.c lnp.h serial.c serial.h ${TEST} ${BENCH} ${INCLUDES} gaussian_c.o gaussian_s.o gaussian_p.o gaussian_b.o
	${CPP} ${CFLAGS} commit.c lnp.c serial.c shuffle.c sha224-256.c gaussian_c.o gaussian_s.o gaussian_p.o gaussian_b.o ${RAND} ${TEST} ${BENCH} -o $@ ${LIBS}

lnp: lnp.c lnp.h commit.c commit.h ${TEST} ${BENCH} ${INCLUDES} gaussian_c.o gaussian_s.o gaussian_p.o gaussian_b.o
	${CPP} ${CFLAGS} -DLNP_MAIN commit.c lnp.c sha224-256.c gaussian_c.o gaussian_s.o gaussian_p.o gaussian_b.o ${RAND} ${TEST} ${BENCH} -o $@ ${LIBS}

clean:
	rm -f *.o commit encrypt vericrypt shuffle lnp
