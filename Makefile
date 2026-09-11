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
DYNBIN := bin/simphy-x-dynamic

# default: bin/simphy-x is the fully static x86-64 Linux binary (glibc/libm, GSL, MPFR, GMP, SQLite
# baked in): runs on any x86-64 Linux without installing anything and gives identical output on
# every machine.  `make dynamic` builds bin/simphy-x-dynamic against the host's glibc instead.
all: $(BIN)

$(BIN): $(OBJ)
	@mkdir -p bin
	$(CC) $(CFLAGS) $(LDFLAGS) -static $^ -o $@ -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -lgsl -lgslcblas -lsqlite3 -lmpfr -lgmp -lm -ldl 2>&1 | grep -v "dlopen' in statically linked\|unixDlOpen" || true
	@test -x $@ && echo "built $@ (static, portable)"

static: $(BIN)

dynamic: $(DYNBIN)

$(DYNBIN): $(OBJ)
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

.PHONY: all clean tests check install static dynamic
