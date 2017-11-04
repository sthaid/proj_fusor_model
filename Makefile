TARGETS = model

# XXX -Wno-unused-variable

CC = gcc
OUTPUT_OPTION=-MMD -MP -o $@
CFLAGS = -g -O2 -pthread -mcmodel=medium -Wall -Iutil $(shell sdl2-config --cflags) 

SRC_MODEL = main.c \
            model.c \
            util/util_sdl.c \
            util/util_sdl_predefined_panes.c \
            util/util_jpeg.c \
            util/util_png.c \
            util/util_misc.c
OBJ_MODEL=$(SRC_MODEL:.c=.o)

DEP=$(SRC_MODEL:.c=.d)

#
# build rules
#

all: $(TARGETS)

model: $(OBJ_MODEL) 
	$(CC) -pthread -lrt -lm -lpng -ljpeg -lSDL2 -lSDL2_ttf -lSDL2_mixer \
              -o $@ $(OBJ_MODEL)

-include $(DEP)

#
# clean rule
#

clean:
	rm -f $(TARGETS) $(OBJ_MODEL)$(DEP)

