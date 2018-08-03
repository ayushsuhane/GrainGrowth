IDIR = ../include

#CC=gcc -W -Wall -ggdb3 -g
CC=gcc -O3
CFLAGS=-I$(IDIR)

ODIR = obj
LDIR = ../lib

LIBS= -lm

_DEPS = Global.h Functions.h nodefunctions.h analysisfunctions.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

 
_OBJ = main.o initialize.o node.o pfsolver.o output.o analysis.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) 

GG: $(OBJ) 
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~

