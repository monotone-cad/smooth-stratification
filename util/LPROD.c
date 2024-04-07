#include "../smooth_stratification.h"

// list product
Word LPROD(Word L)
{
    Word m = 1, a = 0;
    while (L != NIL) {
        ADV(L, &a, &L);

        // multiply by zero short-circuits
        if (a == 0) return 0;

        // otherwise perform the multiplication
        m *= a;
    }

    return m;
}


