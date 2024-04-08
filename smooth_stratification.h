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
Word ISEMPTY(Word r, Word L, Word P, Word V);

// input/output
int read_input(Word *r, Word *V, Word *P);
void write_output(Word r, Word S, Word V);
void write_polynomials(Word r, Word Ps, Word V);

// algorithm
Word JacobiFromMinor(Word r, Word P, Word j, Word Hs, Word Is, Word Minor);
Word strat_helper(Word k, Word np, Word r, Word Fs, Word Is, Word Hs, Word Minor, Word *S_, Word V);
Word stratify(Word r, Word Ps, Word *S_, Word V);

#endif
