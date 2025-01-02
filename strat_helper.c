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
void construct_stratum_basic(Word k, Word r, Word V, Word Hs, Word Q, Word Qs, Word Gs, Word Ineqs, Word *k1_, Word *Y_)
{
    Word L, P, i;
    Word Eqs = NIL, Ineqats = NIL;

    // we may have Q = const, if so don't include it.
    if (!IPCONST(r, Q)) {
        Qs = COMP(Q, Qs);
    }

    // add Qs and Hs to the definition of stratum
    L = Hs, i = k - 1;
    while (L != NIL) {
        // equation h = 0
        ADV(L, &P, &L);
        Eqs = COMP(P, Eqs);

       --i;
    }

    L = Qs;
    while (L != NIL) {
        // inequation s /= 0
        ADV(L, &P, &L);
        Ineqats = COMP(P, Ineqats);
    }

    // we start with functions (h_1,...,h_{k-1}) and then determine which ones from Gs are required.
    Word k1 = k - 1; // store the codimension.

    // attempt to add mroe polynomials from the list Gs of candidate functions.
    while (Gs != NIL) {
        ADV(Gs, &P, &Gs);

        if (!ISEMPTY(r, V, Hs, COMP(P, Qs), Ineqs)) {
            // if this set is non-empty, then there is a point at which P /= 0, thus P = 0 is necessary
            Hs = COMP(P, Hs);
            Eqs = COMP(P, Eqs);

            ++k1; // addition of a new polynomial increases the codimension
        }
    }

    // assign to return
    *k1_ = k1;
    *Y_ = LIST2(Eqs, Ineqats);
}

// recursive smooth stratification of polynomials Fs
// np: number of polynomials
// r: number of variables
// Fs: list of polynomials and degrees (..., P_i, deg(P_i), ...)
// Is: indices (i1,...,ik)
// Hs: list of polynomials h_i1,...,h_ik
// Minor: (LENGTH(I1) - 1) * (LENGTH(I1) - 1) - matrix |a_ij| = partial h_i / partial x_j for h in Hs and i in Is
// return list of all differentials computed by alrogithm
Word strat_helper(Word r, Word V, Word Ineqs, Word k, Word np, Word Fs, Word Is, Word Hs, Word Qs, Word Minor, int *strat_count_, Word *S_)
{
    // add constraints x_1 = 0, ..., x_{i1 - 1}  0
    Word i0 = FIRST(Is); // number of variables considered so far

    // base case for no polynomials?
    if (np == 0 || i0 >= r) {
#ifdef DEBUG
        printf("base case, no polynomials or no more derivatives possible.\n");
#endif
        *strat_count_ = 0;
        return NIL;
    }

#ifdef DEBUG
    printf("\nrecursive call, k = %d, i0 = %d\n", k, i0);
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
            if (d == 0) continue;
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
        if (P == 0) break;

        // construct jacobi matrix using h1 = P and i1 = v
        Word Jacobi = JacobiFromMinor(r, P, v, Hs, Is, Minor);

        // compute partial differential, determinant of jacobi matrix
        // note that Q may be a constant, we need to preserve it
        Word Q = MAIPDE(r, Jacobi); // next derivative is the jacobi determinant
        Word Qdeg = DEG(r,Q);
        bool q_const = IPCONST(r, Q);
        Word Q1 = LIST2(Q, Qdeg);
        Word Qs1 = Qs;
        if (!q_const) Qs1 = COMP(Q, Qs1); // because c /= 0 is trivially true

#ifdef DEBUG
        SWRITE("h_k = ");  IPWRITE(r, P, V); SWRITE("\n");
        SWRITE("s_k = ");  IPWRITE(r, Q, V); SWRITE("\n");
#endif

        // candidate stratum Y1 on which Q /= 0
        //   - 0 /= 0 is trivially false, so we immediately conclude that it's empty
        //   - const /= 0 is trivially true, so assuming the algebraic set Gs3 is non-empty, Y1 is trivially non-empty
        //   - if Q is constant, Y1 must be the last candidate, since the next one includes a trivially false equation
        if (Q != 0 && (q_const || !ISEMPTY(r, V, Gs3, Qs1, Ineqs))) {
            // Gs2 contains derivatives computed during recursion
            int strata_appended;
            Gs2 = strat_helper(r, V, Ineqs, k + 1, g_count, Gs, COMP(v, Is), COMP(P, Hs), Qs1, Jacobi, &strata_appended, S_);
            Gs1 = CONC(Gs1, Gs2);

            // determine if all derivatives at step k+1 vanish on the set Y1
            strat_count += strata_appended;

            // append stratum if none were appended during induction
            if (strata_appended == 0) {
                Word Y, k1;
                construct_stratum_basic(k, r, V, Hs1, Q, Qs, Gs3, Ineqs, &k1, &Y);
#ifdef DEBUG
                printf("appending stratum, k = %d, think it has codimension %d\n", k, k1);
#endif
                append_stratum(S_, k1, Y);
                ++strat_count;
            }

            Gs3 = COMP(Q, Gs3);
            Gs = COMP(Q1, Gs);
            ++g_count;

            if (q_const) {
                // the next candidate will include the trivially false eequation const = 0. stop.
                break;
            }
        }

        // append derivative to Fs, preserving zeroes
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

