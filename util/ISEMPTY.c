#include "../smooth_stratification.h"

/**
 * is_empty(r, L, V, Q)
 *
 * cad-based emptiness check - based on SYSSOLVECAD
 * r : number of variables
 * L : list of polynomials
 * V : list of variables
 * return : true if empty, false otherwise
 */
Word ISEMPTY(Word r, Word L, Word V)
{
    Word F = NIL, P, Ct, Cf;
    while (L != NIL) {
        ADV(L, &P, &L);

        F = COMP(LIST4(EQOP, FIRST(P), r, NIL), F);
    }

    F = COMP(ANDOP, F);
    LWRITE(F); SWRITE("\n");

    // re-initialise qepcad before each run
    QepcadCls Q;
	INITSYS();

    // set input formula
    Q.SETINPUTFORMULA(V,LIST4(r, r, NIL, F));
    // Q.PRDQFF();
    Q.CADautoConst();

    // special case: trivially false
    if (Q.GVPC == 0) {
        return 0;
    }

    LISTOFCWTV(Q.GVPC, &Ct, &Cf);

    /* compute cad */
    return Ct == NIL;
}

