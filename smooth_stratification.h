#ifndef __SMOOTH_STRATIFICATION__
#define __SMOOTH_STRATIFICATION__

#include <stdio.h>
#include <stdlib.h>
#include "replacesac.h"
#include "constants.h"
#include "qepcad.h"

// utility functions
Word DEG(Word r, Word P);
Word LSUM(Word L);
Word LPROD(Word L);
Word ISEMPTY(Word r, Word L, Word P, Word Qs, Word V);

// input/output
int read_input(Word *r, Word *V, Word *P, Word *Q);
void write_polynomial(Word r, Word P, Word V);
void write_output(Word r, Word S, Word Qs, Word V);
void write_polynomials(Word r, Word Ps, Word V);

// algorithm
Word JacobiFromMinor(Word r, Word P, Word j, Word Hs, Word Is, Word Minor);
Word strat_helper(Word r, Word V, Word Ineqs, Word k, Word np, Word Fs, Word Is, Word Hs, Word Minor, int *strat_count, Word *S_);
Word stratify(Word r, Word Ps, Word Ineqs, Word V, Word *S_);

#endif
