#ifndef __SMOOTH_STRATIFICATION__
#define __SMOOTH_STRATIFICATION__

// #include <readline/readline.h>
// #include <readline/history.h>

#include "replacesac.h"
#include "constants.h"

// utility functions
Word DEG(Word r, Word P);
Word LSUM(Word L);
Word LPROD(Word L);

// input/output
int read_input(Word *r, Word *V, Word *P);
void write_output(Word r, Word V, Word *S);
void write_polynomials(Word r, Word Ps, Word V);

// algorithm
Word JacobiFromMinor(Word r, Word P, Word j, Word Hs, Word Is, Word Minor);
Word strat_helper(Word k, Word np, Word r, Word Fs, Word Is, Word Hs, Word Minor, Word V);
Word stratify(Word r, Word Ps, Word V);

#endif
