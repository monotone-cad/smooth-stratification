#include "smooth_stratification.h"

void write_polynomials(Word r, Word Ps, Word V)
{
    while (Ps != NIL) {
        Word P;
        ADV(Ps, &P, &Ps);

        write_polynomial(r, P, V);
        SWRITE("\n");
    }
}

