
SRC_C   = $(wildcard *.c)
DIR_C   = $(notdir $(SRC_C))
OBJ_C  = $(patsubst %.c,%.o,$(DIR_C))
BIN_C  = $(patsubst %.c,%,$(DIR_C))

default: lib $(BIN_C)

lib:
	@(cd ../; make)

include ../Makefile.inc

arms: arms.o
iluc: iluc.o
iluk: iluk.o
ilut: ilut.o
vbiluk: vbiluk.o
vbilut: vbilut.o
solver: solver.o

clean:
	rm -f *.o $(BIN_C)
