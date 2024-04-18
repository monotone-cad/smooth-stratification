#include "../smooth_stratification.h"

/**
 * is_empty(r, L, P1, V, Q)
 *
 * cad-based emptiness check - based on SYSSOLVECAD
 * r : number of variables
 * Fs: list of polynomial equations
 * Gs: list of polynomials inequations
 * Hs: list of atomic QEPCAD formulas for strict polynomial inequalities g > 0
 * V : list of variables
 * return : true if empty, false otherwise
 */
Word ISEMPTY(Word r, Word Fs, Word Gs, Word Ineqs, Word V)
{
    Word P, Ct, Cf;

    // start with the list of inequalities
    Word F = Ineqs;

    // add equations
    while (Fs != NIL) {
        ADV(Fs, &P, &Fs);

        F = COMP(LIST4(EQOP, P, r, NIL), F);
    }

    // and inequations
    while (Gs != NIL) {
        ADV(Gs, &P, &Gs);

        F = COMP(LIST4(NEOP, P, r, NIL), F);
    }

    // complete formula by adding the inequation and conujnction
    F = COMP(ANDOP, F);

    // re-initialise qepcad before each run
    QepcadCls Q;
	INITSYS();

    // set input formula
    Q.SETINPUTFORMULA(V,LIST4(r, r, NIL, F));
    //Q.PRDQFF();
    Q.CADautoConst();

    // special case: trivially false
    if (Q.GVPC == 0) {
        return true;
    }

    LISTOFCWTV(Q.GVPC, &Ct, &Cf);

    /* compute cad */
    return Ct == NIL;
}

