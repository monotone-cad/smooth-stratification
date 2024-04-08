#include "smooth_stratification.h"

// perform smooth stratification
Word stratify(Word r, Word L, Word *S_, Word V)
{
    Word D, P;

    Word Fs = NIL, s = 0;
    while (L != NIL) {
        ADV(L, &P, &L);
        ++s;

        D = DEG(r, P);
        Fs = COMP(LIST2(P, D), Fs);
    }

    // initial i0 = FIRST(I1) = 0. h0 = FIRST(Hs) = 0, Minor is the empty matrix
    int strata_appended;
    Word Gs = strat_helper(1, s, r, Fs, LIST1(0), LIST1(0), NIL, &strata_appended, S_, V);

    if (strata_appended == 0) {
        // X is smooth
        printf("TODO X is smooth!\n");
    }

    return Gs;
}

