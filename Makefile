TARGETS = t1 t2 

CC = gcc
OUTPUT_OPTION=-MMD -MP -o $@
CFLAGS = -c -g -O2 -pthread -fsigned-char -Wall \
         $(shell sdl2-config --cflags) 

SRC_T1 = t1.c 
OBJ_T1=$(SRC_T1:.c=.o)

SRC_T2 = t2.c \
         util_misc.c
OBJ_T2=$(SRC_T2:.c=.o)

DEP=$(SRC_T2:.c=.d) $(SRC_T1:.c=.d)

#
# build rules
#

all: $(TARGETS)

t1: $(OBJ_T1) 
	$(CC) -o $@ $(OBJ_T1)

# XXX pg
t2: $(OBJ_T2) 
	$(CC) -pg -lm -o $@ $(OBJ_T2)

-include $(DEP)

#
# clean rule
#

clean:
	rm -f $(TARGETS) $(OBJ_T1) $(OBJ_T2) $(DEP)

