#include <stdio.h>
#include <stdlib.h>
#include "smooth_stratification.h"

int main(int argc, char **argv)
{
    /* Start up SACLIB. read data from stdin, robustly, continuing until read is successful. */
    Word r, Ps, V, *S, t = 0;
    Word ac; char **av;
    ARGSACLIB(argc,argv,&ac,&av);
    BEGINSACLIB((Word *)&argc);

    while (t != 1) {
        t = read_input(&r, &V, &Ps);
    }

    /* Write out the polynomial for testing purposes */
    SWRITE("DEBUG: polynomials entered:\n");
    Word Ps1 = Ps, P;
    while (Ps1 != NIL) {
        ADV(Ps1, &P, &Ps1);
        IPWRITE(r, P, V);
        SWRITE("\n");
    }
    SWRITE("\n");

    stratify(r, Ps, V);

    /* Write output */
    // write_output(r, V, S);

    /* tidy up */
    ENDSACLIB(SAC_FREEMEM);
    // free(S); /* strata array */
    free(av); /* saclib arguments */

    return 0;
}

