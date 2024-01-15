#ifndef __SMOOTH_STRATIFICATION__
#define __SMOOTH_STRATIFICATION__

// #include <readline/readline.h>
// #include <readline/history.h>

#include "replacesac.h"
#include "constants.h"

int read_input(Word *r, Word *V, Word *P);

Word stratify(Word r, Word Ps);

void write_output(Word r, Word V, Word *S);

#endif
