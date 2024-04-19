#include "smooth_stratification.h"

// perform smooth stratification
Word stratify(Word r, Word L, Word Ineqs, Word V, Word *S_)
{
    // initialise strata S
    *S_ = NIL;
    int i = 0;
    while (i < r) {
        *S_ = COMP(NIL, *S_);
        ++i;
    }

    // initialise the input set of polynomials with their degrees
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
    Word Gs = strat_helper(r, V, Ineqs, 1, s, Fs, LIST1(0), LIST1(0), NIL, NIL, &strata_appended, S_);

    if (strata_appended == 0) {
        // X is smooth
        *S_ = NIL;
    }

    return Gs;
}

