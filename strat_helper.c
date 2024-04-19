#include "smooth_stratification.h"

// convenience function:
// append a stratum Y, codimension k.
void append_stratum(Word *S_, Word k, Word Y)
{
    Word S = LELTI(*S_, k);
    S = COMP(Y, S);
    SLELTI(*S_, k, S);
}

// construct a stratum in basic representation.
void construct_stratum_basic(Word k, Word r, Word V, Word Hs, Word Q, Word Gs, Word Ineqs, Word *k1_, Word *Y_)
{
    Word P;
    Word Qs = IPCONST(r,Q) ? NIL : LIST1(Q); // inequations

    // candidate list of polynomials begins with Hs. we then see if any of Gs need to be appended.
    --k; // Hs so far contains h_1,...,h_{k-1}
    Word k1 = k; // codimension >= k
    while (Gs != NIL) {
        ADV(Gs, &P, &Gs);

        if (!ISEMPTY(r, V, Hs, COMP(P, Qs), Ineqs)) {
            // if this set is non-empty, then there is a point at which P /= 0, thus P adds some information.
            Hs = COMP(P, Hs);
            ++k1;
        }
    }

    // assign k1
    *k1_ = k1;

    // construct the stratum and determine the codimension
    // begin with Q /= 0
    Word Y = LIST3(LIST2('s', k+1), Q, NEOP);

    // add Hs: extra polynomials first
    while (k1 > k && Hs != NIL) {
        ADV(Hs, &P, &Hs);
        Y = COMP3(NIL, P, EQOP, Y);
        --k1;
    }

    // add the remaining Hs
    Word i = 0;
    while (i < k && Hs != NIL) {
        ++i;
        ADV(Hs, &P, &Hs);
        Y = COMP3(LIST2('h', i), P, EQOP, Y);
    }

    *Y_ = Y;
}

// recursive smooth stratification of polynomials Fs
// np: number of polynomials
// r: number of variables
// Fs: list of polynomials and degrees (..., P_i, deg(P_i), ...)
// Is: indices (i1,...,ik)
// Hs: list of polynomials h_i1,...,h_ik
// Minor: (LENGTH(I1) - 1) * (LENGTH(I1) - 1) - matrix |a_ij| = partial h_i / partial x_j for h in Hs and i in Is
// return list of all differentials computed by alrogithm
Word strat_helper(Word r, Word V, Word Ineqs, Word k, Word np, Word Fs, Word Is, Word Hs, Word Minor, int *strat_count_, Word *S_)
{
    // add constraints x_1 = 0, ..., x_{i1 - 1}  0
    Word i0 = FIRST(Is); // number of variables considered so far
    Word v;

    // polynomial 1 x_v^1, in r variables.
    if (i0 > 0) {
        Ineqs = COMP(LIST4(EQOP, PMONSV(r, 1, i0, 1), r, NIL), Ineqs);
    }

    // base case for no polynomials?
    if (np == 0 || k > r) {
#ifdef DEBUG
        printf("base case, no polynomials or k exceeds number of variables. nothing to do\n");
#endif
        *strat_count_ = 0;
        return NIL;
    }

#ifdef DEBUG
    printf("\nrecursive call, k = %d\n", k);
#endif

    // list of H polynomials (rev order) without the last (zero) one.
    Word Hs1 = RED(CINV(Hs));

    // set up return value
    Word Gs1 = NIL; // list of all differentials computed in this round, to return
    Word Gs2 = NIL; // list of differentials produced during induction
    Word Gs3 = NIL; // store Gs on current round

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
        Gs3 = COMP(FIRST(F1), Gs3);

        // increment index
        ++p_index;
    }

    // base case, i1 = r.
    if (i0 == r && !ISEMPTY(r, V, Gs3, NIL, Ineqs)) {
        Word k1, Y;
        construct_stratum_basic(k, r, V, Hs1, 0, Gs3, Ineqs, &k1, &Y);
#ifdef DEBUG
        printf("base case, ik = r but stratum non-empty.\n");
        printf("appending stratum, k = %d, think it has codimension %d\n", k, k1);
#endif
        append_stratum(S_, k, Y);

        *strat_count_ = 1;
        return NIL;
    } else if (i0 == r) { // base case with empty stratum
#ifdef DEBUG
        printf("base case, ik = r but stratum is empty.\n");
#endif

        *strat_count_ = 0;
        return NIL;

    }

    // main loop, consider each index (j, m_{i0 + 1}, ..., m_r) in lex order
    Word n_finished = 0; // each polynomial has a different max index, keep track of how many indices are maxed out
    Word count = 0; // number of derivatives computed so far, scalar value of (m_{i0 + 1}, ..., m_r)

    int strat_count = 0;
    while (n_finished < np) { // stop once differential index for every polynomial is maxed
        // reached the end of polynomial list, cycle back to beginning and consider next differential index I
        if (p_index == np) {
            p_index = 0;
            n_finished = 0;

            ++count;
        }

        v = Dvs[p_index]; // differentiation variable
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

        // candidate stratum Y1 is non-empty
        if (Q != 0 && !ISEMPTY(r, V, Gs3, LIST1(Q), Ineqs)) {
            // Gs2 contains derivatives computed during recursion
            int strata_appended;
            Gs2 = strat_helper(r, V, Ineqs, k + 1, g_count, Gs, COMP(v, Is), COMP(P, Hs), Jacobi, &strata_appended, S_);
            Gs1 = CONC(Gs1, Gs2);

            // determine if all derivatives at step k+1 vanish on the set Y1
            strat_count += strata_appended;

            // append stratum if none were appended during induction
            if (strata_appended == 0) {
                Word Y, k1;
                construct_stratum_basic(k, r, V, Hs1, Q, Gs3, Ineqs, &k1, &Y);
#ifdef DEBUG
                printf("appending stratum, k = %d, think it has codimension %d\n", k, k1);
#endif
                append_stratum(S_, k, Y);
                ++strat_count;
            }

            Gs = COMP(LIST2(Q, Qdeg), Gs);
            Gs3 = COMP(Q, Gs3);
            ++g_count;
        }

        // append derivative to Fs, preserving zeroes
        Word Q1 = LIST2(Q, Qdeg);
        SRED(F1, Q1);
        Append[p_index] = RED2(F1);

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

    Gs1 = CONC(Gs1, Gs3);
    *strat_count_ = strat_count;
    return Gs1;
}

