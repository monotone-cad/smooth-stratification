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

// transform QEPCAD formula into polynomial list.
Word process_formula(Word F, Word r)
{
    Word junk, op, P, F1, Ps = NIL;
    ADV4(F, &junk, &junk, &junk, &F, &junk);
    ADV(F, &op, &F);

    if (op == EQOP) { // single polynomial
        return LIST1(FIRST(F));
    } else if (op != ANDOP) { // only accept conjunctions
        return NIL;
    }

    while (F != NIL) {
        ADV(F, &F1, &F);
        FIRST2(F1, &op, &P);

        if (op != EQOP) {
            return 0;
        }

        Ps = COMP(P, Ps);
    }

    return INV(Ps);
}

int read_input(Word *r_, Word *V_, Word *Ps_)
{
    Word t, r, V, F, Ps = NIL;

	/* Read in variable list for the polynomial */
	SWRITE("Enter a variable list.\n");

	VLREADR(&V, &t);
    if (t == 0) {
        SWRITE("ERROR: invalid variable list.\n");

        return t;
    }

	r = LENGTH(V);

    /* Read in polynomials as a QEEPCAD formula. */
    SWRITE("Please enter a QEPCAD formula defining an elementary semialgebraic set");
    SWRITE("(i.e., conjunction of equations and inequalities).\n");

    FREADR(V, r, &F, &t);
    if (t == 0) {
        SWRITE("ERROR: invalid QEPCAD formula.\n");

        return t;
    }

    // transform formula into a list of polynomials
    Ps = process_formula(F, r);
    LWRITE(Ps); SWRITE("\n");
    if (Ps == NIL) {
        SWRITE("ERROR: invalid formula. Please enter a conjunction of polynomial equations.");

        return 0;
    }

    *r_ = r;
    *V_ = V;
    *Ps_ = Ps;
    return 1;
}

