#include "smooth_stratification.h"

/*======================================================================
 * Write out computed strata
 *      r: number of variables
 *      V: variable list, [x1,...,xr]
 *      S: list of strata
 *  Ineqs: list of inequalities
 *
 *====================================================================*/

void write_op(Word op)
{
    switch (op) {
        case LTOP:
            SWRITE("<  0");
            break;
        case EQOP:
            SWRITE("== 0");
            break;
        case GTOP:
            SWRITE(">  0");
            break;
        case GEOP:
            SWRITE(">= 0");
            break;
        case NEOP:
            SWRITE("/= 0");
            break;
        case LEOP:
            SWRITE("<= 0");
            break;
        default:
            SWRITE("???");
    }
}

void write_label(Word L)
{
    // no label
    if (L == NIL) {
        SWRITE("   ");

        return;
    }

    Word pref, ind;
    FIRST2(L, &pref, &ind);

    CWRITE(pref);
    SWRITE("_");
    IWRITE(ind);
}

void write_output(Word r, Word S, Word Ineqs, Word V)
{
    // special case: input set is smooth
    if (S == NIL) {
        SWRITE("Stratification algorithm determined that input set X is smooth.\n");

        return;
    }

    if (Ineqs != NIL) {
        SWRITE("Strata are contained in an open set defined by the following strict inequalities.\n");
    }

    while (Ineqs != NIL) {
        Word P;
        ADV(Ineqs, &P, &Ineqs);

        // Note: P is an atomic QEPCAD formula
        SWRITE("  "); write_polynomial(r, SECOND(P), V); SWRITE(" "); write_op(GTOP); SWRITE("\n");
    }

    Word k = 0;
    while (S != NIL) {
        Word Sk, Sk1;
        ADV(S, &Sk, &S);
        ++k;

        SWRITE("Strata of codimension "); GWRITE(k); SWRITE("\n");

        // strata, codimension k.
        Word j = 0;
        while (Sk != NIL) {
            ADV(Sk, &Sk1, &Sk);
            ++j;

            // write label
            SWRITE("X_("); IWRITE(k); SWRITE(","); IWRITE(j); SWRITE("):\n");

            // Write stratum
            SWRITE("  { ");
            while (Sk1 != NIL) {
                Word L, P, op;
                ADV3(Sk1, &L, &P, &op, &Sk1);

                write_label(L); SWRITE(" ");
                write_polynomial(r, P, V); SWRITE(" "); write_op(op); SWRITE("\n  ");

                if (Sk1 == NIL) {
                    SWRITE("} ");
                } else {
                    SWRITE(", ");
                }
            }
            SWRITE("\n");
        }

        SWRITE("\n");
    }
}

