# simphy-x build
CC ?= gcc
DEPS ?= $(CURDIR)/deps
CFLAGS ?= -O2 -std=gnu11 -Wall -Wextra -Wno-unused-parameter -ffp-contract=off -fno-fast-math
CPPFLAGS := -I$(DEPS)/usr/include -I$(DEPS)/usr/include/x86_64-linux-gnu
LDFLAGS := -L$(DEPS)/usr/lib/x86_64-linux-gnu
STATIC_LIBS := -Wl,-Bstatic -lgsl -lgslcblas -lsqlite3 -lmpfr -lgmp -Wl,-Bdynamic
LIBS := $(STATIC_LIBS) -lm -ldl -lpthread

SRC := $(wildcard src/*.c)
OBJ := $(patsubst src/%.c,build/%.o,$(SRC))
BIN := bin/simphy-x

all: $(BIN)

# fully static, portable x86-64 Linux binary (libc/libm baked in: identical output on any machine)
static: $(OBJ)
	@mkdir -p bin
	$(CC) $(CFLAGS) $(LDFLAGS) -static $^ -o bin/simphy-x-static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -lgsl -lgslcblas -lsqlite3 -lmpfr -lgmp -lm -ldl
	@echo "built bin/simphy-x-static"

$(BIN): $(OBJ)
	@mkdir -p bin
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LIBS)

build/%.o: src/%.c src/*.h
	@mkdir -p build
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

TESTSRC := src/sx_unif.c src/sx_coal.c src/sx_common.c src/sx_rng.c
TESTLIBS := -Wl,-Bstatic -lgsl -lgslcblas -lmpfr -lgmp -Wl,-Bdynamic -lm -lpthread

build/unif_test: tests/unif_test.c $(TESTSRC) src/*.h
	@mkdir -p build
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $< $(TESTSRC) -o $@ $(TESTLIBS)
build/bridge_exact_test: tests/bridge_exact_test.c $(TESTSRC) src/*.h
	@mkdir -p build
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $< $(TESTSRC) -o $@ $(TESTLIBS)
build/bridge_vs_inversion: tests/bridge_vs_inversion.c $(TESTSRC) src/*.h
	@mkdir -p build
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $< $(TESTSRC) -o $@ $(TESTLIBS)
build/simphy_cdf_check: tests/simphy_cdf_check.c $(TESTSRC) src/*.h
	@mkdir -p build
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $< $(TESTSRC) -o $@ $(TESTLIBS)
build/compat_prob_test: tests/compat_prob_test.c $(TESTSRC) src/*.h
	@mkdir -p build
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $< $(TESTSRC) -o $@ $(TESTLIBS)

tests: build/unif_test build/bridge_exact_test build/bridge_vs_inversion build/simphy_cdf_check build/compat_prob_test

# numerics unit test + byte-exact matrix against the reference build (needs upstream/SimPhy/bin/simphy_ref)
check: $(BIN) build/unif_test
	./build/unif_test
	tests/verify_compat.sh /tmp/simphy-x-verify 1 4

install: $(BIN)

	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 755 $(BIN) $(DESTDIR)$(PREFIX)/bin/simphy-x

PREFIX ?= /usr/local

clean:
	rm -rf build bin

.PHONY: all clean tests check install static
