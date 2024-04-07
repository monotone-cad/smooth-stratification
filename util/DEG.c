#include "../smooth_stratification.h"

// convenience function for multi-degree of polynomial, returns (d_r,...,d_1) where d_i is degree of P in variable i
// returns degree + 1
Word DEG(Word r, Word P)
{
    Word D = PDEGV(r, P);
    Word d, D1 = NIL;

    while (D != NIL) {
        ADV(D, &d, &D);

        D1 = COMP(d+1, D1);
    }

    return D1;
}


