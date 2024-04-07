#include "smooth_stratification.h"

// construct a jacobi (k + 1) * (k + 1)-matrix from a minor (k*k-matrix) by appending one column
// (\H_k / \x_j,...,\H_1) and one row (\P / \x_j, \P / \x_ik,..., \P / \x_i1)
Word JacobiFromMinor(Word r, Word P, Word j, Word Hs, Word Is, Word Minor)
{
    // Step 1: add column
    Word Jacobi = NIL; // list of k rows
    Word k = 0; // dimension of Minor

    Word H, Row;
    while (Minor != NIL) {
        ADV(Minor, &Row, &Minor);
        ADV(Hs, &H, &Hs);
        ++k;

        // compute \H / \x_j
        Word H1 = IPDER(r, H, j);

        // append H1 to row
        Row = COMP(H1, Row);
        Jacobi = COMP(Row, Jacobi);
    }

    // Jacobi was constructed backwards
    Jacobi = INV(Jacobi);

    Word i;
    Row = NIL;
    while (k > 0) { // Is is assumed to have k+1 elements
        ADV(Is, &i, &Is);
        --k;

        // compute \P / \x_i and append
        Word P1 = IPDER(r, P, i);
        Row = COMP(P1, Row);
    }

    // row constructed in order (i1,...,ik)
    Row = INV(Row);

    // compute \P / \x_j and append
    Word P1 = IPDER(r, P, j);
    Row = COMP(P1, Row);

    // add row to jacobi matrix and return
    return COMP(Row, Jacobi);
}


