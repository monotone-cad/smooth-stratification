#include "smooth_stratification.h"

void write_polynomials(Word r, Word Ps, Word V)
{
    while (Ps != NIL) {
        Word P;
        ADV(Ps, &P, &Ps);

        IPWRITE(r, P, V);
        SWRITE("\n");
    }
}

