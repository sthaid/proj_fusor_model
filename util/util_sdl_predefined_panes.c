/*
Copyright (c) 2017 Steven Haid

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdarg.h>
#include <stddef.h>
#include <unistd.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>
#include <limits.h>
#include <assert.h>

#include "util_sdl.h"
#include "util_sdl_predefined_panes.h"
#include "util_misc.h"

// -----------------  PREDEFINED PANE HANDLERS  -----------------------------

int32_t pane_handler_display_text(pane_cx_t * pane_cx, int32_t request, void * init, sdl_event_t * event) 
{
    #define PIXELS_PER_ROW   (sdl_font_char_height(0))

    #define SDL_EVENT_USER_MOUSE_MOTION   (SDL_EVENT_USER_DEFINED + 0)   // xxx better name then EVENT_USER
    #define SDL_EVENT_USER_MOUSE_WHEEL    (SDL_EVENT_USER_DEFINED + 1)

    struct {
        texture_t texture[1000];
        int32_t max_texture;
        int32_t text_y;
    } * vars = pane_cx->vars;

    // ----------------------------
    // -------- INITIALIZE --------
    // ----------------------------

    if (request == PANE_HANDLER_REQ_INITIALIZE) {
        char *newline, *p;
        char  line[200];
        char *text = init;

        vars = pane_cx->vars = calloc(1,sizeof(*vars));
        for (p = text; *p; ) {
            newline = strchr(p, '\n');
            if (newline) {
                memcpy(line, p, newline-p);
                line[newline-p] = '\0';
                p = newline + 1;
            } else {
                strcpy(line, p);
                p += strlen(line);
            }
            vars->texture[vars->max_texture++] = sdl_create_text_texture(WHITE, BLACK, 0, line);
        }
        return PANE_HANDLER_RET_NO_ACTION;
    }

    // ------------------------
    // -------- RENDER --------
    // ------------------------

    if (request == PANE_HANDLER_REQ_RENDER) {
        rect_t * pane = &pane_cx->pane;
        int32_t  lines = pane->h / PIXELS_PER_ROW;
        int32_t  i;

        // sanitize text_y, this is the location of the text that is displayed
        // at the top of the display
        if (vars->text_y > PIXELS_PER_ROW * (vars->max_texture - lines + 2)) {
            vars->text_y = PIXELS_PER_ROW * (vars->max_texture - lines + 2);
        }
        if (vars->text_y < 0) {
            vars->text_y = 0;
        } 

        // display the text
        for (i = 0; i < vars->max_texture; i++) {
            sdl_render_texture(pane, 0, i*PIXELS_PER_ROW-vars->text_y, vars->texture[i]);
        }

        // register control events 
        rect_t loc = {0,0,pane->w,pane->h};
        sdl_register_event(pane, &loc, SDL_EVENT_USER_MOUSE_MOTION, SDL_EVENT_TYPE_MOUSE_MOTION, pane_cx);
        sdl_register_event(pane, &loc, SDL_EVENT_USER_MOUSE_WHEEL, SDL_EVENT_TYPE_MOUSE_WHEEL, pane_cx);

        return PANE_HANDLER_RET_NO_ACTION;
    }

    // -----------------------
    // -------- EVENT --------
    // -----------------------

    if (request == PANE_HANDLER_REQ_EVENT) {
        int32_t  lines = pane_cx->pane.h / PIXELS_PER_ROW;

        switch (event->event_id) {
        case SDL_EVENT_USER_MOUSE_MOTION:
            vars->text_y -= event->mouse_motion.delta_y;
            return PANE_HANDLER_RET_DISPLAY_REDRAW;
        case SDL_EVENT_USER_MOUSE_WHEEL:
            vars->text_y -= event->mouse_wheel.delta_y * 2 * PIXELS_PER_ROW;
            return PANE_HANDLER_RET_DISPLAY_REDRAW;
        case SDL_EVENT_KEY_HOME:
            vars->text_y = 0;
            return PANE_HANDLER_RET_DISPLAY_REDRAW;
        case SDL_EVENT_KEY_END:
            vars->text_y = INT_MAX;
            return PANE_HANDLER_RET_DISPLAY_REDRAW;
        case SDL_EVENT_KEY_PGUP:
            vars->text_y -= (lines - 2) * PIXELS_PER_ROW;
            return PANE_HANDLER_RET_DISPLAY_REDRAW;
        case SDL_EVENT_KEY_PGDN:
            vars->text_y += (lines - 2) * PIXELS_PER_ROW;
            return PANE_HANDLER_RET_DISPLAY_REDRAW;
        case SDL_EVENT_KEY_UP_ARROW:
            vars->text_y -= PIXELS_PER_ROW;
            return PANE_HANDLER_RET_DISPLAY_REDRAW;
        case SDL_EVENT_KEY_DOWN_ARROW:
            vars->text_y += PIXELS_PER_ROW;
            return PANE_HANDLER_RET_DISPLAY_REDRAW;
        }

        return PANE_HANDLER_RET_NO_ACTION;
    }

    // ---------------------------
    // -------- TERMINATE --------
    // ---------------------------

    if (request == PANE_HANDLER_REQ_TERMINATE) {
        int32_t i;
        for (i = 0; i < vars->max_texture; i++) {
            sdl_destroy_texture(vars->texture[i]);
        }
        free(vars);
        return PANE_HANDLER_RET_NO_ACTION;
    }

    // not reached
    assert(0);
    return PANE_HANDLER_RET_NO_ACTION;
}

