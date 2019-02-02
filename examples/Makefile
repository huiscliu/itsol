
SRC_C   = $(wildcard *.c)
DIR_C   = $(notdir $(SRC_C))
OBJ_C  = $(patsubst %.c,%,$(DIR_C))
BIN_C  = $(patsubst %.c,%,$(DIR_C))

default: lib $(OBJ_C) $(BIN_C)

lib:
	@(cd ../; make)

include ../Makefile.inc


clean:
	rm -f *.o
