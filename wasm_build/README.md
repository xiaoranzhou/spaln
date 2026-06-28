# spaln.wasm — Emscripten WebAssembly build

`spaln.wasm` + `spaln.js` are a self-contained WebAssembly port of SPALN 3.0.8
compiled with Emscripten 4.0.20.  `spaln.js` is the Emscripten loader/glue;
`spaln.wasm` is the compiled binary.  Both are gitignored (built artifacts);
only the `Makefile` is tracked.

---

## Build

Prerequisites: emsdk activated (`source /path/to/emsdk/emsdk_env.sh`).

```bash
cd spaln-source/wasm_build
make          # produces spaln.js (63 KB) + spaln.wasm (542 KB)
make test     # smoke-test: checks "SPALN version" in output
make clean    # remove all *.o *.a spaln.js spaln.wasm
```

### Exact flags used

**Compile flags** (every `.cc` → `.o`):
```
-O2 -std=c++11 -ffp-contract=off
-DNO_SIMD          # exclude SSE4.1/NEON path — see below
-DSINGLE_THREAD    # disable pthreads-dependent code paths
-I../src
-sUSE_ZLIB=1       # Emscripten-bundled zlib
-Wno-deprecated-declarations
```

**Link flags** (`spaln.js` link step only):
```
-sNO_DISABLE_EXCEPTION_CATCHING   # keep C++ exceptions alive in WASM
-sALLOW_MEMORY_GROWTH             # heap grows as needed
-sINITIAL_MEMORY=67108864         # 64 MB initial heap
-sENVIRONMENT=node                # Node.js glue only
-sEXPORTED_FUNCTIONS='["_main"]'
-sEXPORT_ES6=0                    # CommonJS require()-compatible
```

---

## Why `-DNO_SIMD`

`src/fwd2s1_simd.cc` contains the hot alignment loop implemented with SSE4.1
intrinsics (x86) and NEON (ARM64 via sse2neon).  WebAssembly SIMD uses a
different ISA (`wasm_simd128.h`); the existing intrinsic code does not compile
under `emcc`.

The entire body of `fwd2s1_simd.cc` is wrapped in `#ifndef NO_SIMD … #endif`.
Passing `-DNO_SIMD` makes the file compile to an empty translation unit, and
the scalar fallback in `fwd2h1.cc` takes over.  Output is identical to the
SIMD path; throughput is roughly 3–5× lower on long proteins.

A WASM SIMD port (replacing `_mm_*` with `wasm_*`) is the natural Phase 2
performance target.

---

## Memory64 (`spaln_m64.wasm`)

Standard WASM uses 32-bit linear memory (4 GB address space).  SPALN's internal
structures use `long` and pointer-width arithmetic that assumes 64-bit width on
the host; with a 32-bit WASM heap and a large genome index (≥ ~100 MB `.bkp`
file), pointer arithmetic silently truncates, producing wrong alignments or
crashes.

The memory64 build adds `-sMEMORY64=1` at link time, enabling 64-bit WASM
addressing (WASM `memory64` proposal, supported in Node ≥ 20 and Chrome ≥ 119).
This makes pointer-width types truly 64-bit inside WASM, matching the native
ABI.

| Binary | Size | Use case |
|---|---|---|
| `spaln.wasm` | 542 KB | development, smoke tests, small genomes |
| `spaln_m64.wasm` | 550 KB | full-size genome indices (e.g. 514 MB `.bkp`) |

Use `spaln_m64.wasm` when the `.bkp` genome index exceeds ~100 MB.

---

## ALN_TAB staging (required before any alignment)

SPALN reads scoring matrices and splice-site tables from the directory given by
`-T<dir>` (or the compile-time default).  Inside WASM, the host filesystem is
not visible; all table files must be staged into Emscripten MEMFS before calling
`_main`.

Minimal staging pattern (Node.js):

```js
const SpalnModule = require('./spaln.js');
const fs   = require('fs');
const path = require('path');

const mod = await SpalnModule({ noInitialRun: true });

function mkdirp(p) {
  p.split('/').filter(Boolean).reduce((acc, seg) => {
    const cur = acc + '/' + seg;
    try { mod.FS.mkdir(cur); } catch (_) {}
    return cur;
  }, '');
}
function stage(src, dst) {
  mkdirp(dst);
  for (const e of fs.readdirSync(src)) {
    const s = path.join(src, e), d = dst + '/' + e;
    fs.statSync(s).isDirectory() ? stage(s, d)
                                 : mod.FS.writeFile(d, fs.readFileSync(s));
  }
}

// Stage spaln-source/table/ into MEMFS at /table (must happen before _main)
stage('/path/to/spaln-source/table', '/table');

// Stage genome index files and query FASTA
for (const ext of ['.bkp', '.ent', '.grp', '.idx', '.seq'])
  mod.FS.writeFile('/genome' + ext, fs.readFileSync('/host/genome' + ext));
mod.FS.writeFile('/query.faa', fs.readFileSync('/host/query.faa'));

// Run alignment (pass -T/table so SPALN finds the staged tables)
mod.callMain(['-Q7', '-O4,7', '-T/table', '-d', '/genome', '/query.faa']);

// Read output back
const result = mod.FS.readFile('/out.O7', { encoding: 'utf8' });
```

Always pass `-T/table` explicitly so SPALN finds the staged tables regardless
of the compile-time default.

---

## Smoke test

```bash
# Quick sanity check via Makefile
make test
# Expected: PASS: spaln.wasm runs, outputs usage

# Manual
node spaln.js 2>&1 | grep "SPALN version"
# Expected: *** SPALN version 3.0.8 ***
```

For a full alignment round-trip, see `/tmp/spaln_wasm/test_align.js` — it
stages `Due-2_med.*` index files + `query.faa` and calls `mod._main(argc, argv)`.
