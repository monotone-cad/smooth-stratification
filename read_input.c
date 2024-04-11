#include "smooth_stratification.h"

/*======================================================================
 * Input processing
 * t <- read_input(r, V, P)
 *      r: number of variables
 *      V: variable list, [x1,...,xr]
 *      P: integral polynomial in [x1,...,xr]
 *  Ineqs: list of atomic formulas for strict polynomial inequalities g > 0
 *      t: 1 if successful, 0 otherwise
 *
 *====================================================================*/

// transform QEPCAD formula into polynomial list.
// return true if able to process, false otherwise
bool process_formula(Word r, Word F, Word *Ps_, Word *Ineqs_)
{
    Word junk, op, P, F1, Ps = NIL, Ineqs = NIL;
    ADV4(F, &junk, &junk, &junk, &F, &junk);
    ADV(F, &op, &F);

    if (op == EQOP) { // single polynomial
        *Ps_ = LIST1(FIRST(F));
        *Ineqs_ = NIL;

        return true;
    } else if (op != ANDOP) { // only accept conjunctions
        return false;
    }

    while (F != NIL) {
        ADV(F, &F1, &F);
        FIRST2(F1, &op, &P);

        if (op == EQOP) { // f = 0
            Ps = COMP(P, Ps);
        } else if (op == GTOP) { // g > 0
            Ineqs = COMP(LIST4(GTOP, P, r, NIL), Ineqs);
        } else if (op == LTOP) { // g < 0 => -g > 0
            Ineqs = COMP(LIST4(GTOP, IPNEG(r, P), r, NIL), Ineqs);
        } else { // other operations not permitted.
            return false;
        }
    }

    *Ps_ = INV(Ps); // polynomials appear in order they are entered.
    *Ineqs_ = Ineqs; // order of inequalities does not matter.
    return true;
}

int read_input(Word *r_, Word *V_, Word *Ps_, Word *Ineqs_)
{
    Word t, r, V, F;

	/* Read in variable list for the polynomial */
	SWRITE("Enter a variable list.\n");

	VLREADR(&V, &t);
    if (t == 0) {
        SWRITE("ERROR: invalid variable list.\n");

        return 0;
    }

	r = LENGTH(V);

    /* Read in polynomials as a QEEPCAD formula. */
    SWRITE("Please enter a QEPCAD formula defining an elementary semialgebraic set");
    SWRITE("(i.e., conjunction of equations and inequalities).\n");

    FREADR(V, r, &F, &t);
    if (t == 0) {
        SWRITE("ERROR: invalid QEPCAD formula.\n");

        return 0;
    }

    // transform formula into a list of polynomials
    if (!process_formula(r, F, Ps_, Ineqs_)) {
        SWRITE("ERROR: invalid formula. Please enter a conjunction of polynomial equations.");

        return 0;
    }

    *r_ = r;
    *V_ = V;
    // note: Ps_ and Ineqs_ already assigned in process_formula
    return 1;
}

