CCOMP = gcc
#CFLAGS = -O4 -DNDEBUG -DEXCESS_TYPE_LONG -DPRINT_STAT -DCHECK_SOLUTION -Wall -lm
CFLAGS = -g -DPRINT_FLOW -DEXCESS_TYPE_LONG -DPRINT_STAT -DCHECK_SOLUTION -Wall -lm

all: coptreem 
coptreem: main.c funcs.c parser.c timer.c
	$(CCOMP) $(CFLAGS) -o coptreem main.c libm.so 
clean: 
	rm -f coptreem *~
