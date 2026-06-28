EMCC = /Users/xr/metagraph/emsdk/upstream/emscripten/emcc
EMAR = /Users/xr/metagraph/emsdk/upstream/emscripten/emar

SRC  = ../src
CXXFLAGS = -O2 -std=c++11 -ffp-contract=off \
           -DNO_SIMD -DSINGLE_THREAD \
           -I$(SRC) \
           -sUSE_ZLIB=1 \
           -Wno-deprecated-declarations

# clib.a members
CLIB_OBJS = clib.o iolib.o mfile.o sets.o supprime.o

# sblib.a members (no fwd2s1_simd — guarded by NO_SIMD)
SLIB_OBJS = aln2.o dbs.o gaps.o codepot.o divseq.o kmers.o gsinfo.o \
            fwd2b1.o fwd2d1.o fwd2h1.o fwd2s1.o bitpat.o eijunc.o \
            seq.o simmtx.o sqpr.o utilseq.o vmf.o wln.o boyer_moore.o

all: spaln.js

# ── libraries ────────────────────────────────────────────────────────────────
clib.a: $(CLIB_OBJS)
	$(EMAR) rc $@ $^

sblib.a: $(SLIB_OBJS)
	$(EMAR) rc $@ $^

# ── spaln entry point ────────────────────────────────────────────────────────
spaln.js: spaln.o blksrc.o sblib.a clib.a
	$(EMCC) $(CXXFLAGS) \
	  -sNO_DISABLE_EXCEPTION_CATCHING \
	  -sALLOW_MEMORY_GROWTH \
	  -sINITIAL_MEMORY=67108864 \
	  -sENVIRONMENT=node \
	  -sEXPORTED_FUNCTIONS='["_main"]' \
	  -sEXPORT_ES6=0 \
	  spaln.o blksrc.o sblib.a clib.a \
	  -o spaln.js

# ── pattern rule for .cc → .o ─────────────────────────────────────────────
%.o: $(SRC)/%.cc
	$(EMCC) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o *.a spaln.js spaln.wasm

.PHONY: all clean

# ── WASM acceptance test ─────────────────────────────────────────────────────
NODE = /Users/xr/metagraph/emsdk/node/22.16.0_64bit/bin/node

test: spaln.js
	$(NODE) spaln.js 2>&1 | grep -q "SPALN version" && echo "PASS: spaln.wasm runs, outputs usage" || (echo "FAIL" && exit 1)
