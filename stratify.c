#include <stdio.h>
#include <stdlib.h>
#include "smooth_stratification.h"

// perform smooth stratification
Word stratify(Word r, Word L)
{
    Word D, P;

    Word Fs = NIL, k = 0;
    while (L != NIL) {
        ADV(L, &P, &L);
        ++k;

        D = DEG(r, P);
        Fs = COMP(LIST2(P, D), Fs);
    }

    // initial i0 = FIRST(I1) = 0. h0 = FIRST(Hs) = 0, Minor is the empty matrix
    Word Gs = strat_helper(k, r, Fs, LIST1(0), LIST1(0), NIL);

    return Gs;
}

