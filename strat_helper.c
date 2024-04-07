#include "smooth_stratification.h"

Word construct_stratum(Word Backup[], Word k, Word np, Word p_index, Word h_index, Word s_index)
{
    Word n_null = 0, i = 0, j = 0, L = NIL;

    // avoid modifying the lists stored in Backup
    Word Fs1[np];
    while (i < np) {
        Fs1[i] = Backup[i];
        ++i;
    }

    // loop over each polynomial in lex order
    // note that lists may have different lengths, n_null counts how many are null
    // exit when all are null.
    i = -1; // first time through the loop, I will be incremented.
    while (n_null < np) {
        // update indices
        if (i == np - 1) {
            n_null = 0;
            i = 0;
            ++j;
        } else {
            ++i;
        }

        if (Fs1[i] == NIL) {
            ++n_null;

            // skip if null
            continue;
        }

        Word P, junk;
        ADV2(Fs1[i], &P, &junk, &Fs1[i]);

        // skip zero polynomials
        if (P == 0) {
            continue;
        }

        Word Label = NIL;
        Word op = EQOP;

        // label polynomials s_k and h_k, identified by their index
        if (i == p_index) {
            if (j == h_index) {
                Label = LIST2('h',k);
            } else if (j == s_index) {
                Label = LIST2('s',k);
                op = NEOP; // polynomial s is always not equal to 0
            }
        }

        // append the polynomial
        L = COMP3(Label, P, op, L);
    }

    return L;
}

// recursive smooth stratification of polynomials Fs
// np: number of polynomials
// r: number of variables
// Fs: list of polynomials and degrees (..., P_i, deg(P_i), ...)
// Is: indices (i1,...,ik)
// Hs: list of polynomials h_i1,...,h_ik
// Minor: (LENGTH(I1) - 1) * (LENGTH(I1) - 1) - matrix |a_ij| = partial h_i / partial x_j for h in Hs and i in Is
// return list of all differentials computed by alrogithm
Word strat_helper(Word k, Word np, Word r, Word Fs, Word Is, Word Hs, Word Minor, Word *S_, Word V)
{
    // end of recursion
    Word i0 = FIRST(Is); // number of variables considered so far
    if (np == 0 || i0 == r) {
#ifdef DEBUG
        printf("base case: nothing to do\n");
#endif
        return NIL;
    }

#ifdef DEBUG
    printf("\nrecursive call, k = %d\n", k);
#endif

    // set up return value
    Word Gs1 = NIL; // list of all differentials computed in this round, to return
    Word S = LELTI(*S_, k);

    // set up working array
    Word g_count = np; // how many differentials computed so far, index in Gs
    Word Gs = NIL; // list of all differentials computed in this round, working set

    // metadata
    Word Ms[np]; // number of steps before maximum differentiation variable (i1) should be incremented
    Word Dvs[np]; // list of current differentiation variable (i1)
    Word Backup[np]; // first polynomial in the list
    Word Chase[np]; // first element is h1
    Word ChaseIndex[np]; // chase array index of polynomial h1
    Word Append[np]; // pointer to last but one element in list, for quick appending.

    Word p_index = 0; // initial polynomial index, ranges over 0 <= p_index < np

    // initialise metadata for each input polynomial
    while (p_index < np) {
        Word F1;
        ADV(Fs, &F1, &Fs);

        Word D = REDI(SECOND(F1), i0);
        Ms[p_index] = COMP(1, LCOPY(D)); // neet do copy as we modify later
        Dvs[p_index] = i0;
        Backup[p_index] = F1;
        Chase[p_index] = F1;
        ChaseIndex[p_index] = 0;
        Append[p_index] = RED(F1);
        Gs = COMP(LCOPY(F1), Gs);

        // increment index
        ++p_index;
    }

    // main loop, consider each index (j, m_{i0 + 1}, ..., m_r) in lex order
    Word n_finished = 0; // each polynomial has a different max index, keep track of how many indices are maxed out
    Word count = 0; // number of derivatives computed so far, scalar value of (m_{i0 + 1}, ..., m_r)

    bool is_empty = false;
    while (!is_empty && n_finished < np) { // stop once differential index for every polynomial is maxed
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
            Chase[p_index] = Backup[p_index];
            ChaseIndex[p_index] = 0;

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

#ifdef DEBUG
        printf(
            "p_index %d, variable %d, count = %d, chase_index = %d, len = %d\n",
            p_index, v, count, ChaseIndex[p_index], LENGTH(Backup[p_index]) / 2);
#endif

        // compute s_k = partial_{(h_1,...,h_{k-1}),(i_1,...,i_{k-1}),v} h_k
        // get h_k and its degree
        Word P;
        ADV(Chase[p_index], &P, &Chase[p_index]);

        // construct jacobi matrix using h1 = P and i1 = v
        Word Jacobi = JacobiFromMinor(r, P, v, Hs, Is, Minor);

        // compute partial differential, determinant of jacobi matrix
        // note that Q may be a constant, we need to preserve it
        Word Q = MAIPDE(r, Jacobi); // next derivative is the jacobi determinant
        Word Qdeg = DEG(r,Q);

#ifdef DEBUG
        SWRITE("h_k = ");  IPWRITE(r, P, V); SWRITE("\n");
        SWRITE("s_k = ");  IPWRITE(r, Q, V); SWRITE("\n");
#endif

        if (Q != 0) {
            // Gs2 contains derivatives computed during recursion
            Word Gs2 = strat_helper(k + 1, g_count, r, Gs, COMP(v, Is), COMP(P, Hs), Jacobi, S_, V);
            Gs1 = CONC(Gs1, Gs2);

            Gs = COMP(LIST2(Q, Qdeg), Gs);
            ++g_count;
        }

        // append derivative to Fs, preserving zeroes
        Word Q1 = LIST2(Q, Qdeg);
        SRED(F1, Q1);
        Append[p_index] = RED2(F1);

        if (Q != 0) {
            // append strata
#ifdef DEBUG
            printf("appending stratum, k = %d\n", k);
#endif
            S = COMP(construct_stratum(Backup, k, np, p_index, ChaseIndex[p_index], count), S);

            is_empty = ISEMPTY(r, Gs, V);
        }

        // next polynomial please.
        Chase[p_index] = RED(Chase[p_index]);
        ChaseIndex[p_index] = ChaseIndex[p_index] + 1;
        ++p_index;
    }

#ifdef DEBUG
    printf("End of round. printing Fs\n");
    for (int i = 0; i < np; i++) {
        printf("%d: ", i+1);
        Word LL = Backup[i];
        while (LL != NIL) {
            Word PP, II;
            ADV2(LL, &PP, &II, &LL);
            if (II == 0) continue;

            LWRITE(II); SWRITE(": ");
            IPWRITE(r, PP, V); SWRITE(", ");
        }
        SWRITE("\n");
    }

    SWRITE("\n");
#endif

    // construct list Gs1 of all functions in Gs
    while (Gs != NIL) {
        Word G1;
        ADV(Gs, &G1, &Gs);

        // skip constant and zero polynomials
        if (LSUM(SECOND(G1)) == r) continue;

        Gs1 = COMP(FIRST(G1), Gs1);
    }

    // save strata, codimension k
    SLELTI(*S_, k, S);

    return Gs1;
}

