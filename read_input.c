#include "smooth_stratification.h"

/*======================================================================
 * Input processing
 * t <- read_input(r, V, P)
 *      r: number of variables
 *      V: variable list, [x1,...,xr]
 *      P: integral polynomial in [x1,...,xr]
 *      t: 1 if successful, 0 otherwise
 *
 *====================================================================*/


int read_input(Word *r_, Word *V_, Word *Ps_)
{
    Word r, V, P, Ps = NIL, t;

	/* Read in variable list for the polynomial */
	SWRITE("Enter a variable list:\n");
	V = VLREAD();
	r = LENGTH(V);

	if (r < 1) {
		return 0;
	}

    /* Read in a single polynomial */
    SWRITE("Enter integral polynomials in ");
    VLWRITE(V);
    SWRITE(" variables:\n");

    do {
        IPEXPREAD(r, V, &P, &t);
        FILINE();

        if (t == 0) {
            SWRITE("Error: input polynomial must be in at least one variable.\n");

            return 0;
        }

        Ps = COMP(P, Ps);
    } while (LKAHEAD() != '.');

    *r_ = r;
    *V_ = V;
    *Ps_ = Ps;
    return 1;
}

