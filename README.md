# SMOOTH STRATIFICATION

An algorithm to compute smooth stratifications of semialgebraic sets, based on the work of Gabrielov and Vorobjov. This
work is novel because, rather than using quantifier elimination like standard stratification algorithms, this method
computes partial derivatives of higher orders to find strata of progressively smaller dimension.

## Compile

```bash
make
```

## Usage

The program takes input interactively, a bit like QEPCAD. First enter the variable list in the form `(x_1,...,x_n)`,
then enter each polynomial terminated by a full stop and a new line. At the end of input, enter a single line containing a full stop to notify the
program that input has been read.

### For example
```
(x,y,z)
z.
x^2 - y^2.
.
```
