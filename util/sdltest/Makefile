CC = gcc

SDL2_CFLAGS = $(shell sdl2-config --cflags)

SRC_SDLTEST = sdltest.c \
              ../util_sdl.c \
              ../util_sdl_predefined_panes.c \
              ../util_jpeg.c \
              ../util_png.c \
              ../util_misc.c

#
# build rules
#

sdltest: $(SRC_SDLTEST)
	$(CC) -g -O2 -pthread -Wall -I.. \
              $(SDL2_CFLAGS) \
              -lrt -lm -lpng -ljpeg -lSDL2 -lSDL2_ttf -lSDL2_mixer \
              -o $@ \
              $(SRC_SDLTEST)

#
# clean rule
#

clean:
	rm -f sdltest

