#include "../smooth_stratification.h"

/**
 * is_empty(r, L, P1, V, Q)
 *
 * cad-based emptiness check - based on SYSSOLVECAD
 * r : number of variables
 * L : list of polynomials
 * P1: polynomial which must not be equal to 0
 * V : list of variables
 * return : true if empty, false otherwise
 */
Word ISEMPTY(Word r, Word L, Word P1, Word V)
{
    // list L contains polynomials of the form ((P_1, I_1), ..., (P_k, I_k)), discard I
    Word F = NIL, P, Ct, Cf;
    while (L != NIL) {
        ADV(L, &P, &L);

        F = COMP(LIST4(EQOP, FIRST(P), r, NIL), F);
    }

    // complete formula by adding the inequation and conujnction
    F = COMP2(ANDOP, LIST4(NEOP, P1, r, NIL), F);

    // re-initialise qepcad before each run
    QepcadCls Q;
	INITSYS();

    // set input formula
    Q.SETINPUTFORMULA(V,LIST4(r, r, NIL, F));
    // Q.PRDQFF();
    Q.CADautoConst();

    // special case: trivially false
    if (Q.GVPC == 0) {
        return true;
    }

    LISTOFCWTV(Q.GVPC, &Ct, &Cf);

    /* compute cad */
    return Ct == NIL;
}

