# SPALN WASM Build

## Overview

This directory contains the Emscripten/WebAssembly build of SPALN for the browser pipeline. The WASM binary is used by `example4/spaln_cnv_pipeline/index.html` (Pyodide) so that SPALN alignment and genome-index building can run entirely client-side.

Build outputs (tracked by `Makefile`, binaries gitignored):

- `spaln.js` — Emscripten loader/glue
- `spaln.wasm` — compiled SPALN binary

## Build command

```bash
source /Users/xr/metagraph/emsdk/emsdk_env.sh
cd /Users/xr/git/spaln/spaln-source/src
mkdir -p /tmp/spaln_wasm

# Compile all object files with NO_SIMD and MEMORY64
emcc -O2 -std=c++11 -ffp-contract=off -DNO_SIMD -sMEMORY64=1 -sUSE_ZLIB=1 -sUSE_PTHREADS=0 \
  -c clib.cc iolib.cc mfile.cc sets.cc supprime.cc \
     aln2.cc dbs.cc gaps.cc codepot.cc divseq.cc kmers.cc gsinfo.cc \
     fwd2b1.cc fwd2d1.cc fwd2h1.cc fwd2s1.cc fwd2s1_simd.cc \
     bitpat.cc eijunc.cc seq.cc simmtx.cc sqpr.cc utilseq.cc vmf.cc wln.cc \
     boyer_moore.cc blksrc.cc

# Link the WASM module
emcc -O2 -std=c++11 -ffp-contract=off -DNO_SIMD -sMEMORY64=1 -sUSE_ZLIB=1 -sUSE_PTHREADS=0 \
  -sALLOW_MEMORY_GROWTH=1 -sFORCE_FILESYSTEM=1 -sMODULARIZE=1 -sEXPORT_NAME=SpalnModule \
  -sEXPORTED_RUNTIME_METHODS='["FS","UTF8ToString","ccall","callMain"]' \
  -sEXPORTED_FUNCTIONS='["_main","_malloc","_free"]' \
  -sENVIRONMENT=web,node -sNO_EXIT_RUNTIME=1 \
  spaln.cc *.o -lz -o /tmp/spaln_wasm/spaln_m64.js
```

## Why `-DNO_SIMD`

`fwd2s1_simd.cc` assumes SSE4.1 or ARM NEON intrinsics are available. Emscripten does not define `__SSE4_1__` or `__ARM_NEON__` by default, and passing x86-only intrinsics headers fails on macOS ARM hosts. Adding `-DNO_SIMD` disables the SIMD object (the file body is wrapped in `#ifndef NO_SIMD`) and lets the scalar alignment path in `fwd2s1.cc` / `fwd2b1.cc` / `fwd2h1.cc` do the work.

## ALN_TAB / table staging

The Emscripten build does **not** read the host `ALN_TAB` environment variable. SPALN resolves scoring tables relative to its current working directory, which defaults to `/table` in the WASM build.

Before any alignment or index-build call, stage the contents of `spaln-source/table/` into MEMFS at `/table`:

```javascript
mod.FS.mkdirTree('/table');
for (const [relPath, content] of tableFiles) {
    mod.FS.writeFile('/table/' + relPath, content);
}
mod.FS.chdir('/table');   // relative -T <species> now resolves
```

In the Python bridge (`scripts/wasm_spaln.py`) this is handled by `_stage_table_to_memfs()`.

## Smoke test

```bash
node spaln.js
```

Expected: SPALN usage banner and `*** SPALN version ... ***` printed. A full alignment smoke test is:

```bash
node -e "
const SpalnModule = require('./spaln.js');
SpalnModule({noInitialRun:true}).then(async mod => {
  // stage table/, genome index, query as needed
  mod.FS.chdir('/table');
  mod.callMain(['-po','-Q7','-A1','-T','cladgray','-O4,7','-o','/tmp/out','-d','/tmp/index','/tmp/query.faa']);
  console.log('O4 size:', mod.FS.readFile('/tmp/out.O4').length);
});
"
```

Success criterion: `.O4` file is non-empty and contains the SPALN exon-summary header.
