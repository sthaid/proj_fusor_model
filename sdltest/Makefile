TARGETS = sdltest

CC = gcc
OUTPUT_OPTION=-MMD -MP -o $@
CFLAGS = -g -O2 -pthread -Wall -I../util $(shell sdl2-config --cflags) 

SRC_SDLTEST = sdltest.c \
              ../util/util_sdl.c \
              ../util/util_sdl_predefined_panes.c \
              ../util/util_jpeg.c \
              ../util/util_png.c \
              ../util/util_misc.c
OBJ_SDLTEST=$(SRC_SDLTEST:.c=.o)

DEP=$(SRC_SDLTEST:.c=.d)

#
# build rules
#

all: $(TARGETS)

sdltest: $(OBJ_SDLTEST) 
	$(CC) -pthread -lrt -lm -lpng -ljpeg -lSDL2 -lSDL2_ttf -lSDL2_mixer \
              -o $@ $(OBJ_SDLTEST)

-include $(DEP)

#
# clean rule
#

clean:
	rm -f $(TARGETS) $(OBJ_SDLTEST) $(DEP)

