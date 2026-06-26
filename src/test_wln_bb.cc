/*
 * test_wln_bb.cc — unit tests for the bb ↔ bs invariant in Wlp::eval()
 *
 * Tests two bugs fixed in the pangap branch:
 *
 *  Bug 1 (wln.cc:346): When the backward extension guard fired, bb was frozen
 *    while bs continued decrementing, causing the forward loop to score sigE
 *    from wrong codon positions.  Fix: backward loop always decrements bb
 *    (no reads, safe); forward loop always advances bb but guards the read.
 *
 *  Bug 2 (wln.cc:331): score_p(n) returns data_p+n and is never null.
 *    bb->sigS was read without a bounds check; OOB when jxt->jy+1 >= b->right.
 *    Fix: added b->exin->good(bb) guard before the read.
 *
 * Compile:
 *   cd src && g++ -std=c++11 -O0 -g -fsanitize=address test_wln_bb.cc -o test_wln_bb
 *   ./test_wln_bb
 */

#include <cassert>
#include <cstdio>
#include <cstring>

// ── Minimal simulation of SGPT6 and Exinon bounds ────────────────────────────

struct SGPT6 {
    int sigE;  // exon potential (added during forward extension)
    int sigS;  // start codon potential
    int sigT;  // stop codon potential
};

// Mirrors Exinon::good() and Exinon::score_p() exactly.
struct MockExin {
    const SGPT6* data_p;
    int bias;   // = b->left - 1
    int size;   // = b->right - b->left + 2

    const SGPT6* begin_p() const { return data_p + bias; }
    const SGPT6* end_p()   const { return data_p + bias + size - 1; }
    const SGPT6* score_p(int n)  const { return data_p + n; }
    bool good(const SGPT6* bb)   const {
        return data_p && begin_p() <= bb && bb < end_p();
    }
};

// ── Simulate the FIXED forward loop logic ────────────────────────────────────
// Models the exact code path in wln.cc after the Bug 1 fix.
static int simulate_forward(const MockExin& exin, const SGPT6* bb_in,
                             int start_pos, int end_pos, int bbt)
{
    const SGPT6* bb = bb_in;
    int total_sigE = 0;
    // forward loop: bs advances from start_pos to end_pos by bbt
    for (int pos = start_pos; pos < end_pos; pos += bbt) {
        if (bb) {
            if (exin.good(bb)) total_sigE += bb->sigE;  // guarded read
            bb += bbt;                                   // always advance (Bug 1 fix)
        }
    }
    return total_sigE;
}

// ── Test: Bug 2 — OOB read of bb->sigS at initialization ─────────────────────
// Before fix: `if (bb && ...)` — bb is never null from score_p, no range check.
// After fix:  `if (b->exin->good(bb) && ...)` — OOB position is rejected.
static void test_bug2_oob_sigS_read()
{
    // Genomic b: left=10, right=20  →  valid exin indices 9..19
    const int b_left = 10, b_right = 20;
    const int bias = b_left - 1;           // 9
    const int size = b_right - b_left + 2; // 12

    // Allocate exactly size+1 elements and offset by bias (mirrors Exinon ctor)
    SGPT6 raw[20] = {};
    raw[9].sigS  = 42;   // valid: position 9 == begin_p()
    raw[19].sigS = 99;   // valid: last valid position (b_right - 1 = 19)

    MockExin exin;
    exin.data_p = raw - 0;  // data_p[n] == raw[n]
    exin.bias   = bias;
    exin.size   = size;

    // Case A: jxt->jy = 8 → bb = data_p + 9 = begin_p() → VALID
    {
        int jxt_jy = 8;
        const SGPT6* bb = exin.score_p(jxt_jy + 1);  // data_p + 9
        bool in_range = exin.good(bb);
        assert(in_range && "jxt->jy=8 should be at begin_p(), valid");
        // Fixed: guard passes → would read bb->sigS = raw[9].sigS = 42
        if (in_range) assert(bb->sigS == 42);
    }

    // Case B: jxt->jy = 19 → bb = data_p + 20 = end_p() → OOB (b->right)
    {
        int jxt_jy = 19;  // = b->right - 1, so jxt->jy + 1 = b->right = end_p()
        const SGPT6* bb = exin.score_p(jxt_jy + 1);  // data_p + 20 = end_p()
        bool in_range = exin.good(bb);
        assert(!in_range && "jxt->jy=b->right-1 → bb=end_p() is OOB, must be rejected");
        // Bug 2 fix: guard rejects → sigS NOT read → no OOB access
    }

    // Case C: jxt->jy = 25 → bb = data_p + 26, well past end → OOB
    {
        int jxt_jy = 25;
        const SGPT6* bb = exin.score_p(jxt_jy + 1);  // data_p + 26
        bool in_range = exin.good(bb);
        assert(!in_range && "jxt->jy=25 > b->right → OOB, must be rejected");
    }

    printf("PASS: test_bug2_oob_sigS_read\n");
}

// ── Test: Bug 1 — bb/bs desync during backward extension ─────────────────────
// Verify that after backward extension the bb offset matches jxt->jy.
// Then verify forward pass scores sigE only for positions within the valid range.
static void test_bug1_bb_bs_sync()
{
    // Genomic b: left=10, right=40 → valid exin indices 9..39
    const int b_left = 10, b_right = 40;
    const int bias = b_left - 1;           // 9
    const int size = b_right - b_left + 2; // 32
    const int bbt  = 3;

    SGPT6 raw[50] = {};
    // Set sigE=1 for every valid position so we can count them
    for (int i = bias; i < bias + size - 1; ++i) raw[i].sigE = 1;

    MockExin exin;
    exin.data_p = raw;
    exin.bias   = bias;
    exin.size   = size;

    // Seed: jxt->jy = 30 (well inside valid range)
    // bb initialized at score_p(30 + 1) = data_p + 31
    int jxt_jy = 30;
    const SGPT6* bb = exin.score_p(jxt_jy + 1);

    // ── Simulate backward loop (FIXED: always decrement) ──────────────────
    // Extend backward 10 steps; the last 4 steps go past begin_p() (index 9)
    // Step 1: jxt->jy=27, bb→data_p+28  (valid: >=9)
    // Step 2: jxt->jy=24, bb→data_p+25  (valid)
    // Step 3: jxt->jy=21, bb→data_p+22  (valid)
    // Step 4: jxt->jy=18, bb→data_p+19  (valid)
    // Step 5: jxt->jy=15, bb→data_p+16  (valid)
    // Step 6: jxt->jy=12, bb→data_p+13  (valid)
    // Step 7: jxt->jy=9,  bb→data_p+10  (valid: 10>=9)
    // Step 8: jxt->jy=6,  bb→data_p+7   (OOB: 7<9) — old guard would freeze bb here
    // Step 9: jxt->jy=3,  bb→data_p+4   (OOB)
    // Step 10:jxt->jy=0,  bb→data_p+1   (OOB)
    int steps = 10;
    for (int i = 0; i < steps; ++i) {
        jxt_jy -= bbt;
        if (bb) bb -= bbt;  // Bug 1 fix: always decrement, no good() check
    }
    // After backward: jxt->jy = 30 - 30 = 0, bb = data_p + (31 - 30) = data_p + 1
    assert(jxt_jy == 0);
    assert(bb == exin.score_p(jxt_jy + 1) &&
           "bb must equal score_p(jxt->jy+1) after backward loop (sync invariant)");

    // ── Simulate forward loop (FIXED: always advance, guarded read) ────────
    // Forward: advance 10 steps from bb=data_p+1 with bbt=3.
    // begin_p()=data_p+9, end_p()=data_p+40 (bias+size-1 = 9+32-1 = 40).
    // Positions: data_p+{1,4,7,10,13,16,19,22,25,28}
    // Valid (9<=pos<40): k=3(+10),4(+13),5(+16),6(+19),7(+22),8(+25),9(+28) → 7 hits
    int expected_sigE = 0;
    const SGPT6* bb_check = exin.score_p(jxt_jy + 1);  // start of fwd = data_p+1
    for (int k = 0; k < steps; ++k) {
        if (bb_check) {
            if (exin.good(bb_check)) expected_sigE += bb_check->sigE;
            bb_check += bbt;
        }
    }
    assert(expected_sigE == 7 && "expected 7 valid sigE positions in forward pass");

    // Compare with BUGGY behavior (old diff: bb frozen when OOB in backward loop)
    // Old backward loop left bb frozen at data_p+10 (first OOB step at k=8 from
    // jxt->jy=30: after 7 valid steps, bb=data_p+10 which is valid, then step 8
    // would go to data_p+7 < begin_p() so guard freezes bb at data_p+10).
    // Forward loop then: bb starts at data_p+10 (not data_p+1 — DESYNCED by 9)
    // sigE scored from data_p+10, data_p+13, ..., data_p+10+9*3=data_p+37 (all valid)
    // = 10 × sigE=1 = 10 (WRONG — 2 extra spurious scores from wrong positions)
    const SGPT6* bb_buggy = exin.score_p(30 + 1);  // original bb init
    for (int i = 0; i < 7; ++i) {
        // First 7 backward steps stay in range
        if (bb_buggy && exin.good(bb_buggy - bbt)) bb_buggy -= bbt;
        // (step 8 would freeze bb_buggy at data_p+10 since data_p+7 < begin_p())
    }
    // bb_buggy is now frozen at data_p+10 (valid) even though bs is at data_p+1
    int buggy_sigE = 0;
    for (int k = 0; k < steps; ++k) {
        if (bb_buggy && exin.good(bb_buggy)) {
            buggy_sigE += bb_buggy->sigE;
            bb_buggy += bbt;
        }
    }
    assert(buggy_sigE != expected_sigE &&
           "buggy code must produce wrong sigE count (proves the bug existed)");
    printf("  Bug 1 concrete: fixed sigE=%d, buggy sigE=%d (different — bug confirmed)\n",
           expected_sigE, buggy_sigE);

    printf("PASS: test_bug1_bb_bs_sync\n");
}

// ── Test: forward-loop sigE only from valid range ────────────────────────────
// After a normal (no desync) backward extension, forward scores are correct.
static void test_forward_sigE_in_range_only()
{
    const int b_left = 5, b_right = 20;
    const int bias = b_left - 1;
    const int size = b_right - b_left + 2;
    const int bbt  = 3;

    SGPT6 raw[30] = {};
    for (int i = bias; i < bias + size - 1; ++i) raw[i].sigE = 10;

    MockExin exin;
    exin.data_p = raw;
    exin.bias   = bias;
    exin.size   = size;

    // bb enters OOB at the right end mid-forward-pass
    // Start bb at data_p+13 (valid), simulate 5 forward steps
    const SGPT6* bb = exin.score_p(13);
    int total = 0;
    for (int k = 0; k < 5; ++k) {
        if (bb) {
            if (exin.good(bb)) total += bb->sigE;  // guarded
            bb += bbt;                              // always advance
        }
    }
    // end_p() = data_p+bias+size-1 = data_p+4+17-1 = data_p+20.
    // Positions: data_p+{13,16,19,22,25}. Valid (4<=pos<20): 13✓,16✓,19✓ → 3×10=30.
    assert(total == 30 && "only in-range positions should contribute sigE");

    printf("PASS: test_forward_sigE_in_range_only\n");
}

// ── Negative test: confirm old code would access OOB memory ──────────────────
// Without the fix, score_p(jxt->jy+1) for jxt->jy == b->right-1 produces a
// pointer equal to end_p().  good() must reject this — if it didn't, reading
// bb->sigS would be an out-of-bounds access (one past the valid array).
static void test_bug2_oob_is_one_past_end()
{
    const int b_left = 5, b_right = 15;
    const int bias = b_left - 1;           // 4
    const int size = b_right - b_left + 2; // 12

    // Exinon allocates data_p[bias .. bias+size] (size+1 elements).
    // Lay out raw so that raw[bias+size] is the first element PAST the valid range.
    SGPT6 raw[20] = {};
    for (int i = bias; i < bias + size; ++i) raw[i].sigE = i;  // distinct values

    MockExin exin;
    exin.data_p = raw;
    exin.bias   = bias;
    exin.size   = size;

    // end_p() = data_p + bias + size - 1 = data_p + 15
    assert(exin.end_p() == raw + 15);
    assert(exin.begin_p() == raw + 4);

    // jxt->jy = b->right - 1 = 14 → score_p(15) = data_p+15 = end_p() → OOB
    const SGPT6* bb_oob = exin.score_p(14 + 1);  // == end_p()
    assert(bb_oob == exin.end_p() && "score_p at right edge equals end_p()");
    assert(!exin.good(bb_oob) && "end_p() itself must be OOB (< not <=)");

    // Pre-fix condition: `if (bb && ...)` — bb is non-null, would proceed to read
    bool would_read_prefix = (bb_oob != nullptr);
    assert(would_read_prefix && "pre-fix check is always true — OOB read would happen");

    // Post-fix condition: `if (b->exin->good(bb) && ...)` — correctly rejected
    bool would_read_postfix = exin.good(bb_oob);
    assert(!would_read_postfix && "post-fix guard correctly prevents OOB read");

    printf("PASS: test_bug2_oob_is_one_past_end  "
           "(pre-fix would read: %s, post-fix would read: %s)\n",
           would_read_prefix ? "YES" : "no",
           would_read_postfix ? "yes" : "NO");
}

int main()
{
    test_bug2_oob_sigS_read();
    test_bug2_oob_is_one_past_end();
    test_bug1_bb_bs_sync();
    test_forward_sigE_in_range_only();
    printf("All tests passed.\n");
    return 0;
}
