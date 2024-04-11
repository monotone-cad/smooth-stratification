CC=g++
ARFLAGS = rvU

INCLUDES=-I$(saclib)/include -I$(qe)/source

EXTLIBS=-lreadline -lcurses

EXTLIBSOPT=$(qe)/source/qepcad.a \
${qe}/extensions/sfext/sfexto.a \
${qe}/extensions/lift2D/lift2Do.a \
${qe}/extensions/newadj/newadjo.a \
${qe}/extensions/adj2d/adj2do.a \
${qe}/extensions/rend/rendo.a \
${saclib}/lib/saclibo.a \
$(EXTLIBS)

EXTLIBSDEB=$(qe)/source/qepcad.a \
${qe}/extensions/sfext/sfexto.a \
${qe}/extensions/lift2D/lift2Do.a \
${qe}/extensions/newadj/newadjo.a \
${qe}/extensions/adj2d/adj2do.a \
${qe}/extensions/rend/rendo.a \
${saclib}/lib/saclibd.a \
$(EXTLIBS)

CFLAGSBOTH=-Wall
CFLAGSDEB=$(CFLAGSBOTH) -g
CFLAGSDEBPRINT=$(CFLAGSDEB) -DDEBUG
CFLAGSOPT=$(CFLAGSBOTH) -O4

LIBOPT=stratify.a
LIBDEB=stratifyd.a
LIBDEBP=stratifydp.a
EXE=main

DEPENDENCIESOPT=\
$(LIBOPT)(util/DEG.o) \
$(LIBOPT)(util/LPROD.o) \
$(LIBOPT)(util/LSUM.o) \
$(LIBOPT)(util/ISEMPTY.o) \
$(LIBOPT)(JacobiFromMinor.o) \
$(LIBOPT)(strat_helper.o) \
$(LIBOPT)(read_input.o) \
$(LIBOPT)(stratify.o) \
$(LIBOPT)(write_output.o) \
$(LIBOPT)(write_polynomial.o) \
$(LIBOPT)(write_polynomials.o) \
$(LIBOPT)(main.o) \

DEPENDENCIESDEB=\
$(LIBDEB)(util/DEG.o) \
$(LIBDEB)(util/LPROD.o) \
$(LIBDEB)(util/LSUM.o) \
$(LIBDEB)(util/ISEMPTY.o) \
$(LIBDEB)(JacobiFromMinor.o) \
$(LIBDEB)(strat_helper.o) \
$(LIBDEB)(read_input.o) \
$(LIBDEB)(stratify.o) \
$(LIBDEB)(write_output.o) \
$(LIBDEB)(write_polynomial.o) \
$(LIBDEB)(write_polynomials.o) \
$(LIBDEB)(main.o) \

DEPENDENCIESDEBP=\
$(LIBDEBP)(util/DEG.o) \
$(LIBDEBP)(util/LPROD.o) \
$(LIBDEBP)(util/LSUM.o) \
$(LIBDEBP)(util/ISEMPTY.o) \
$(LIBDEBP)(JacobiFromMinor.o) \
$(LIBDEBP)(strat_helper.o) \
$(LIBDEBP)(read_input.o) \
$(LIBDEBP)(stratify.o) \
$(LIBDEBP)(write_output.o) \
$(LIBDEBP)(write_polynomial.o) \
$(LIBDEBP)(write_polynomials.o) \
$(LIBDEBP)(main.o) \


all: opt deb

# optimised
opt:override CFLAGS = $(CFLAGSOPT)
opt:$(DEPENDENCIESOPT)
	ranlib $(LIBOPT)
	$(CC) $(CFLAGSOPT) $(INCLUDES) $(LIBOPT) $(EXTLIBSOPT) $(LIBOPT) $(EXTLIBSOPT) -o $(EXE)

# debug
deb:override CFLAGS = $(CFLAGSDEB)
deb:$(DEPENDENCIESDEB)
	ranlib $(LIBDEB)
	$(CC) $(CFLAGSDEB) $(INCLUDES) \
        $(LIBDEB) $(EXTLIBSDEB) \
        $(LIBDEB) $(EXTLIBSDEB) \
        -o $(EXE)

# debug and print
debprint:override CFLAGS = $(CFLAGSDEBPRINT)
debprint:$(DEPENDENCIESDEBP)
	ranlib $(LIBDEBP)
	$(CC) $(CFLAGSDEBPRINT) $(INCLUDES) \
        $(LIBDEBP) $(EXTLIBSDEB) \
        $(LIBDEBP) $(EXTLIBSDEB) \
        -o $(EXE)


# clean
clean:
	rm -f main $(LIBDEB) $(LIBDEBP) $(LIBOPT)

# objects
# %.o:override CC=gcc
%.o:%.c
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

