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

int32_t pane_hndlr_display_text(pane_cx_t * pane_cx, int32_t request, void * init_params, sdl_event_t * event) 
{
    #define PIXELS_PER_ROW   (sdl_font_char_height(0))

    #define SDL_EVENT_MOUSE_MOTION   (SDL_EVENT_USER_DEFINED + 0)
    #define SDL_EVENT_MOUSE_WHEEL    (SDL_EVENT_USER_DEFINED + 1)

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
        char *text = init_params;

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
        sdl_register_event(pane, &loc, SDL_EVENT_MOUSE_MOTION, SDL_EVENT_TYPE_MOUSE_MOTION, pane_cx);
        sdl_register_event(pane, &loc, SDL_EVENT_MOUSE_WHEEL, SDL_EVENT_TYPE_MOUSE_WHEEL, pane_cx);

        return PANE_HANDLER_RET_NO_ACTION;
    }

    // -----------------------
    // -------- EVENT --------
    // -----------------------

    if (request == PANE_HANDLER_REQ_EVENT) {
        int32_t  lines = pane_cx->pane.h / PIXELS_PER_ROW;

        switch (event->event_id) {
        case SDL_EVENT_MOUSE_MOTION:
            vars->text_y -= event->mouse_motion.delta_y;
            return PANE_HANDLER_RET_DISPLAY_REDRAW;
        case SDL_EVENT_MOUSE_WHEEL:
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

// XXX move this to util_sdl
// XXX need smaller font
// XXX review and finish up
int32_t pane_hndlr_display_graph(pane_cx_t * pane_cx, int32_t request, void * init_params, sdl_event_t * event) 
{
    #define SDL_EVENT_MOUSE_MOTION   (SDL_EVENT_USER_DEFINED + 0)
    #define SDL_EVENT_MOUSE_WHEEL    (SDL_EVENT_USER_DEFINED + 1)
    #define SDL_EVENT_RESET_GRAPH    (SDL_EVENT_USER_DEFINED + 2)

    struct {
        pane_hndlr_display_graph_params_t ** params;
        float y_min_orig;
        float y_max_orig;
        float y_min;
        float y_max;
    } * vars = pane_cx->vars;
    rect_t * pane = &pane_cx->pane;

    // ----------------------------
    // -------- INITIALIZE --------
    // ----------------------------

    if (request == PANE_HANDLER_REQ_INITIALIZE) {
        vars = pane_cx->vars = calloc(1,sizeof(*vars));
        vars->params    = (pane_hndlr_display_graph_params_t**)init_params;
        return PANE_HANDLER_RET_NO_ACTION;
    }

    // ------------------------
    // -------- RENDER --------
    // ------------------------

    if (request == PANE_HANDLER_REQ_RENDER) {
        int32_t xo, xm, yo, ym, title_col, i;  // xxx variable names
        pane_hndlr_display_graph_params_t * params = *vars->params;

        // if no params then return
        if (params == NULL) {
            return PANE_HANDLER_RET_NO_ACTION;
        }

        // if params y_min/max don't match vars y_min/max then XXX comment
        // reset vars y_min/max (this should not usually happen)
        if (params->y_min != vars->y_min_orig || params->y_max != vars->y_max_orig) {
            INFO("XXXXXXXXXXXXXXXXXX CALLER RESET ymin/max  %f %f\n", params->y_min, params->y_max);
            vars->y_min_orig = params->y_min;
            vars->y_max_orig = params->y_max;
            vars->y_min = params->y_min;
            vars->y_max = params->y_max;
// XXX treet xmin max same
        }

        // init variables that define the location of the x,y axis
        xo = 6 * sdl_font_char_width(0) + 2;  // xxx names
        xm = pane->w - 1;
        yo = pane->h - sdl_font_char_height(0) - 2;
        ym = sdl_font_char_height(0);

        // render the x/y axis, title, units, and scale
//locate y intercept for x == 0
//only if in range
        if (vars->y_min <= 0 && vars->y_max >= 0) {
            int32_t YYY =  yo + (0 - vars->y_min) * ((ym - yo) / (vars->y_max - vars->y_min));
            sdl_render_line(pane, xo, YYY, xm, YYY, WHITE);
        }

        sdl_render_line(pane, xo, yo, xo, ym, WHITE);

        title_col = (pane->w / sdl_font_char_width(0) - strlen(params->title_str)) / 2;
        sdl_render_printf( pane, 0, title_col, 0, WHITE, BLACK, "%s", params->title_str);
        // XXX etc

        // render the graph
        point_t points[params->max_points];
        int32_t max_points = 0;
        for (i = 0; i < params->max_points; i++) {  // XXX vvv optimize
            float x = params->points[i].x;
            float y = params->points[i].y;
            if (x < params->x_min || x > params->x_max ||
                y < vars->y_min || y > vars->y_max) 
            {
                sdl_render_lines(pane, points, max_points, WHITE);
                max_points = 0;
                continue;
            }
            points[max_points].x =  xo + (x - params->x_min) * ((xm - xo) / (params->x_max - params->x_min));
            points[max_points].y =  yo + (y - vars->y_min) * ((ym - yo) / (vars->y_max - vars->y_min));
            max_points++;
        }
        sdl_render_lines(pane, points, max_points, WHITE);

        // register control events 
        rect_t loc = {0,0,pane->w,pane->h};
        sdl_register_event(pane, &loc, SDL_EVENT_MOUSE_MOTION, SDL_EVENT_TYPE_MOUSE_MOTION, pane_cx);
        sdl_register_event(pane, &loc, SDL_EVENT_MOUSE_WHEEL, SDL_EVENT_TYPE_MOUSE_WHEEL, pane_cx);
        sdl_render_text_and_register_event(pane, 0, -1, 0, "R", LIGHT_BLUE, BLACK,
                                           SDL_EVENT_RESET_GRAPH, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);


        return PANE_HANDLER_RET_NO_ACTION;
    }

    // -----------------------
    // -------- EVENT --------
    // -----------------------

    // xxx change y scale with wheel
    // xxx and a RESET button to reset the scale
    if (request == PANE_HANDLER_REQ_EVENT) {
        switch (event->event_id) {
        case SDL_EVENT_MOUSE_MOTION: {
            int32_t yo = pane->h - sdl_font_char_height(0) - 2;
            int32_t ym = sdl_font_char_height(0);
            float delta = -(vars->y_max - vars->y_min) / (ym - yo); 
            vars->y_min += (event->mouse_motion.delta_y * delta);
            vars->y_max += (event->mouse_motion.delta_y * delta);
            return PANE_HANDLER_RET_DISPLAY_REDRAW; }
        case SDL_EVENT_MOUSE_WHEEL: {
            float factor = (event->mouse_motion.delta_y > 0 ? 1.1 : (1. / 1.1));
            //if (vars->y_min <= 0 && vars->y_max >= 0) {
                vars->y_min *= factor;
                vars->y_max *= factor;
            //} else {
                //float avg = (vars->y_min + vars->y_max) / 2.;
                //vars->y_max += (vars->y_max - avg) * factor;
                //vars->y_min += (vars->y_min - avg) * factor;
            //}
            return PANE_HANDLER_RET_DISPLAY_REDRAW; }
        case SDL_EVENT_RESET_GRAPH: {
            vars->y_min = vars->y_min_orig;
            vars->y_max = vars->y_max_orig;
            return PANE_HANDLER_RET_DISPLAY_REDRAW; }
// XXX need RESET
        }
        return PANE_HANDLER_RET_NO_ACTION;
    }

    // ---------------------------
    // -------- TERMINATE --------
    // ---------------------------

    if (request == PANE_HANDLER_REQ_TERMINATE) {
        free(vars);
        return PANE_HANDLER_RET_NO_ACTION;
    }

    // not reached
    assert(0);
    return PANE_HANDLER_RET_NO_ACTION;
}
