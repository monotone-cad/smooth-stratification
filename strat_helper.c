#include <stdio.h>
#include <stdlib.h>
#include "smooth_stratification.h"

// recover index from count
// TODO this is only for debugging, can be deleted when done.
// count = (c1 + c2 * m1 + ... + cn * m(n-1))
Word IndexHelper(Word m, Word M, Word* count_)
{
    Word count = *count_;
    Word a, M1;
    ADV(M, &a, &M1);

    if (M1 == NIL) {
        *count_ = count % m;

        return LIST1(count / m);
    }

    // skip round
    if (a == 1) {
        return COMP(0, IndexHelper(m, M1, count_));
    }

    Word I1 = IndexHelper(m * a, M1, count_);
    count = *count_;

    Word j = count / m;
    *count_ = count % m;

    return COMP(j, I1);
}

// calls index helper
inline Word INDEX(Word count, Word M)
{
    if (count < 0) return NIL;

    Word J = IndexHelper(1, M, &count);

    return J;
}

// recursive smooth stratification of polynomials Fs
// np: number of polynomials
// r: number of variables
// Fs: list of polynomials and degrees (..., P_i, deg(P_i), ...)
// Is: indices (i1,...,ik)
// Hs: list of polynomials h_i1,...,h_ik
// Minor: (LENGTH(I1) - 1) * (LENGTH(I1) - 1) - matrix |a_ij| = partial h_i / partial x_j for h in Hs and i in Is
// return list of all differentials computed by alrogithm
Word strat_helper(Word np, Word r, Word Fs, Word Is, Word Hs, Word Minor)
{
    // end of recursion
    Word i0 = FIRST(Is); // number of variables considered so far
    if (np == 0 || i0 == r) {
        return NIL;
    }

    // TODO debugging
    // SWRITE("strat_helper: Is = "); LWRITE(Is); SWRITE("\n");
    // SWRITE("       Hs = "); LWRITE(Hs); SWRITE("\n");

    // set up return value
    Word Gs1 = NIL; // list of all differentials computed in this round, to return

    // set up working array
    Word g_count = 0; // how many differentials computed so far, index in Gs
    Word Gs = NIL; // list of all differentials computed in this round, working set

    // metadata
    Word Degrees[np]; // degree of Fs[0][0], gives max index
    Word Ms[np]; // number of steps before maximum differentiation variable (i1) should be incremented
    Word Dvs[np]; // list of current differentiation variable (i1)
    Word ChaseDebug[np]; // chase array index of polynomial h1
    Word Backup[np]; // first polynomial in the list
    Word Chase[np]; // first element is h1
    Word Append[np]; // pointer to last but one element in list, for quick appending.

    Word p_index = 0; // initial polynomial index, ranges over 0 <= p_index < np

    // initialise metadata for each input polynomial
    while (p_index < np) {
        Word F1;
        ADV(Fs, &F1, &Fs);

        Word D = REDI(SECOND(F1), i0);

        Degrees[p_index] = D;
        Ms[p_index] = COMP(1, LCOPY(D)); // neet do copy as we modify later
        Dvs[p_index] = i0;
        ChaseDebug[p_index] = 0;
        Backup[p_index] = F1;
        Chase[p_index] = F1;
        Append[p_index] = RED(F1);

        // increment index
        ++p_index;
    }

    // main loop, consider each index (j, m_{i0 + 1}, ..., m_r) in lex order
    Word n_finished = 0; // each polynomial has a different max index, keep track of how many indices are maxed out
    Word count = 0; // number of derivatives computed so far, scalar value of (m_{i0 + 1}, ..., m_r)

    while (n_finished < np) { // stop once differential index for every polynomial is maxed
        // reached the end of polynomial list, cycle back to beginning and consider next differential index I
        if (p_index == np) {
            p_index = 0;
            n_finished = 0;

            ++count;
        }

        Word v = Dvs[p_index]; // differentiation variable
        Word m = FIRST(Ms[p_index]);

        // update variable v and chaser list
        if (count >= m && v == r) { // rollover, but we're finished
            ++n_finished; // this polynomial is done.
            ++p_index; // next polynomial

            continue; // skip it
        } else if (count >= m) { // rollover - increment differentiation variable
            ++v; // next variable ...
            Dvs[p_index] = v; // ... and store
            ChaseDebug[p_index] = 0;
            Chase[p_index] = Backup[p_index];

            // calculate next m
            Word M1 = RED(Ms[p_index]);
            Word d = FIRST(M1);
            SFIRST(M1, d * m);
            Ms[p_index] = M1;

            // degree zero - no derivatives taken for this variable. next iteration will increment the variable.
            if (d == 1) continue;
        }
                // next polynomial
        Word F1 = Append[p_index];

        printf("p_index %d, variable %d, count = %d, chase_index = %d, len = %d\n", p_index, v, count,
        ChaseDebug[p_index], LENGTH(Backup[p_index]));

        // compute s_k = partial_{(h_1,...,h_{k-1}),(i_1,...,i_{k-1}),v} h_k
        // get h_k and its degree
        Word D = Degrees[p_index];
        Word P;
        LWRITE(Chase[p_index]), SWRITE("\n");
        ADV(Chase[p_index], &P, &Chase[p_index]);

        // construct jacobi matrix using h1 = P and i1 = v
        Word Jacobi = JacobiFromMinor(r, P, v, Hs, Is, Minor);

        // compute partial differential, determinant of jacobi matrix
        Word Q = MAIPDE(r, Jacobi); // next derivative is the jacobi determinant
        Word Qdeg = DEG(r,Q);

        if (LSUM(Qdeg) == r) { // s_k is constant, may as well be zero
            Qdeg = 0;
            Q = 0;
        }

        // TODO debugging
        IWRITE(count); SWRITE(", variable: "); IWRITE(v); SWRITE(" polynomial: "); IWRITE(p_index);
        SWRITE(", degree: "); LWRITE(Degrees[p_index]); SWRITE("\n  ");
        LWRITE(INDEX(count, D)); SWRITE(" ");
        P == 0 ? SWRITE("0") : LWRITE(P); SWRITE("\n  ");
        LWRITE(INDEX(ChaseDebug[p_index], D)); SWRITE(" ");
        Q == 0 ? SWRITE("0") : LWRITE(Q); SWRITE("\n\n");

        if (Q != 0) {
            // Gs2 contains derivatives computed during recursion
            Word Gs2 = strat_helper(g_count, r, Gs, COMP(v, Is), COMP(P, Hs), Jacobi);
            Gs1 = CONC(Gs1, Gs2);

            // TODO
            Gs = COMP(LIST2(Q, Qdeg), Gs);
            ++g_count;
        }

        // append derivative to Fs, preserving zeroes
        Word Q1 = LIST2(Q, Qdeg);
        SRED(F1, Q1);
        Append[p_index] = RED2(F1);

        // next polynomial please.
        ChaseDebug[p_index] = ChaseDebug[p_index] + 1;
        Chase[p_index] = RED(Chase[p_index]);
        ++p_index;
    }

    printf("End of round. printing Fs\n");
    for (int i = 0; i < np; i++) {
        printf("%d: ", i+1);
        LWRITE(Backup[i]);
        SWRITE("\n");
    }

    // construct list Gs1 of all functions in Gs
    while (Gs != NIL) {
        Word G1;
        ADV(Gs, &G1, &Gs);

        Gs1 = COMP(FIRST(G1), Gs1);
    }

    return Gs1;
}

