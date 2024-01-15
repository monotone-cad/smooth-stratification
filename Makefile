CC=g++
ARFLAGS = rvU

INCLUDES=-I$(saclib)/include

EXTLIBS=-lreadline

EXTLIBSOPT=${saclib}/lib/saclibo.a \
$(EXTLIBS)

EXTLIBSDEB=${saclib}/lib/saclibd.a \
$(EXTLIBS)

CFLAGSBOTH=-Wall
CFLAGSDEB=$(CFLAGSBOTH) -g
CFLAGSOPT=$(CFLAGSBOTH) -O4

LIBOPT=stratify.a
LIBDEB=stratifyd.a
EXE=main

DEPENDENCIESOPT=\
$(LIBOPT)(read_input.o) \
$(LIBOPT)(stratify.o) \
$(LIBOPT)(write_output.o) \
$(LIBOPT)(main.o) \

DEPENDENCIESDEB=\
$(LIBDEB)(read_input.o) \
$(LIBDEB)(stratify.o) \
$(LIBDEB)(write_output.o) \
$(LIBDEB)(main.o) \


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

# clean
clean:
	rm -f main $(LIBDEB) $(LIBOPT)

# objects
# %.o:override CC=gcc
%.o:%.c
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

