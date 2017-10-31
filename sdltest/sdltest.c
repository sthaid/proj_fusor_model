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

// defines

#define DEFAULT_WIN_WIDTH  1900
#define DEFAULT_WIN_HEIGHT 1000

// prototypes

static int32_t pane_handler_demo(pane_cx_t * cx, int32_t request, void * init, sdl_event_t * event);
static int32_t pane_handler_points_test(pane_cx_t * pane_cx, int32_t request, void * init, sdl_event_t * event);
static bool display_redraw_needed(uint64_t time_render_us);

// -----------------  MAIN  ------------------------------------------------

int main(int argc, char **argv)
{
    int32_t win_width, win_height;

    win_width  = DEFAULT_WIN_WIDTH;
    win_height = DEFAULT_WIN_HEIGHT;
    if (sdl_init(&win_width, &win_height, true, false) < 0) {
        FATAL("sdl_init %dx%d failed\n", win_width, win_height);
    }

    sdl_pane_manager(display_redraw_needed, 1,
                     pane_handler_demo, NULL, 0, 0, 400, 400, PANE_BORDER_STYLE_STANDARD);

    return 0;
}

// -----------------  PANE HANDLERS  ---------------------------------------

static char * text = "\
this is line 1\n\
this is a really looooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooong line\n\
and here is line 3\n\
and here is line 4\n\
and here is line 5\n\
and here is line 6\n\
and here is line 7\n\
and here is line 8\n\
and here is line 9\n\
and here is line 10\n\
and here is line 11\n\
and here is line 12\n\
and here is line 13\n\
and here is line 14\n\
and here is line 15\n\
and here is line 16\n\
and here is line 17\n\
and here is line 18\n\
and here is line 19\n\
and here is line 20\n\
and here is line 21\n\
and here is line 22\n\
and here is line 23\n\
and here is line 24\n\
and here is line 25\n\
";

static int32_t pane_handler_demo(pane_cx_t * pane_cx, int32_t request, void * init, sdl_event_t * event) 
{
    #define SDL_EVENT_MOUSE_CLICK      (SDL_EVENT_USER_DEFINED+0)
    #define SDL_EVENT_MOUSE_WHEEL      (SDL_EVENT_USER_DEFINED+1)
    #define SDL_EVENT_NEW_DEMO_PANE    (SDL_EVENT_USER_DEFINED+2)
    #define SDL_EVENT_NEW_TEXT_PANE    (SDL_EVENT_USER_DEFINED+3)
    #define SDL_EVENT_NEW_POINTS_PANE  (SDL_EVENT_USER_DEFINED+4)
    #define SDL_EVENT_NEW_DISPLAY      (SDL_EVENT_USER_DEFINED+5)

    #define CIRCLE_RADIUS (pane_cx->pane.w / 10)

    struct {
        int32_t   instance;
        uint64_t  time_of_last_render_call_us;
        double    circle_x;
        double    circle_y;
        double    circle_x_rate;
        double    circle_y_rate;
        texture_t circle_texture;
        int32_t   update_counter;
        int32_t   mouse_click_test_counter;
        int32_t   mouse_wheel_test_counter;
        int32_t   key_event_id;
    } * vars = pane_cx->vars;

    static int32_t instance;

    // ----------------------------
    // -------- INITIALIZE --------
    // ----------------------------

    if (request == PANE_HANDLER_REQ_INITIALIZE) {
        vars = pane_cx->vars = calloc(1,sizeof(*vars));
        vars->instance = ++instance;
        vars->time_of_last_render_call_us = 0;
        vars->circle_x = pane_cx->pane.w / 2;
        vars->circle_y = pane_cx->pane.h / 2;
        vars->circle_x_rate = 200;   // pixels per second
        vars->circle_y_rate = 100;   // pixels per second
        vars->circle_texture = sdl_create_filled_circle_texture(CIRCLE_RADIUS, ORANGE);
        vars->update_counter = 0;
        vars->mouse_click_test_counter = 0;
        vars->mouse_wheel_test_counter = 0;
        vars->key_event_id = 0; 
        return PANE_HANDLER_RET_NO_ACTION;
    }

    // ------------------------
    // -------- RENDER --------
    // ------------------------

    if (request == PANE_HANDLER_REQ_RENDER) {
        rect_t * pane = &pane_cx->pane;
        uint64_t time_now_us;
        double   time_delta_secs;

        // get the current time, and init time_of_last_render_call_us
        time_now_us = microsec_timer();
        if (vars->time_of_last_render_call_us == 0) {
            vars->time_of_last_render_call_us = time_now_us;
        }

        // display update_counter and pane instance
        vars->update_counter++;
        sdl_render_printf(pane, 0, 0, 0, WHITE, BLACK, "%d %d", 
            vars->update_counter, vars->instance);

        // keyboard test
        sdl_render_printf(pane, 1, 0, 0, WHITE, BLACK, "KEY 0x%x '%c'", 
            vars->key_event_id, vars->key_event_id ? vars->key_event_id : ' ');

        // 'MOUSE_CLICK' will increment counter on click
        sdl_render_text_and_register_event(
            pane, 3, 0, 0, "MOUSE_CLICK", LIGHT_BLUE, BLACK, 
            SDL_EVENT_MOUSE_CLICK, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);
        sdl_render_printf(pane, 3, 12, 0, WHITE, BLACK, "%d", vars->mouse_click_test_counter);

        // 'MOUSE_WHEEL' will adjust counter on mouse wheel 
        sdl_render_text_and_register_event(
            pane, 4, 0, 0, "MOUSE_WHEEL", LIGHT_BLUE, BLACK, 
            SDL_EVENT_MOUSE_WHEEL, SDL_EVENT_TYPE_MOUSE_WHEEL, pane_cx);
        sdl_render_printf(pane, 4, 12, 0, WHITE, BLACK, "%d", vars->mouse_wheel_test_counter);

        // 'NEW_DEMO_PANE' create another demo pane
        sdl_render_text_and_register_event(
            pane, 5, 0, 0, "NEW_DEMO_PANE", LIGHT_BLUE, BLACK, 
            SDL_EVENT_NEW_DEMO_PANE, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);

        // 'NEW_TEXT_PANE' create another text pane
        sdl_render_text_and_register_event(
            pane, 6, 0, 0, "NEW_TEXT_PANE", LIGHT_BLUE, BLACK, 
            SDL_EVENT_NEW_TEXT_PANE, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);

        // 'NEW_POINTS' create the points test pane
        sdl_render_text_and_register_event(
            pane, 7, 0, 0, "NEW_POINTS", LIGHT_BLUE, BLACK, 
            SDL_EVENT_NEW_POINTS_PANE, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);

        // 'NEW_DISPLAY' chain to another display
        sdl_render_text_and_register_event(
            pane, 8, 0, 0, "NEW_DISPLAY", LIGHT_BLUE, BLACK, 
            SDL_EVENT_NEW_DISPLAY, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);

        // moving circle
        time_delta_secs = (time_now_us - vars->time_of_last_render_call_us) / 1000000.;
        vars->circle_x += vars->circle_x_rate * time_delta_secs;
        if (vars->circle_x < 0) {
            vars->circle_x_rate = -vars->circle_x_rate;
            vars->circle_x = 0;
        }
        if (vars->circle_x > pane->w) {
            vars->circle_x_rate = -vars->circle_x_rate;
            vars->circle_x = pane->w;
        }
        vars->circle_y += vars->circle_y_rate * time_delta_secs;
        if (vars->circle_y < 0) {
            vars->circle_y_rate = -vars->circle_y_rate;
            vars->circle_y = 0;
        }
        if (vars->circle_y > pane->h) {
            vars->circle_y_rate = -vars->circle_y_rate;
            vars->circle_y = pane->h;
        }
        sdl_render_texture(pane, vars->circle_x-CIRCLE_RADIUS, vars->circle_y-CIRCLE_RADIUS, vars->circle_texture);

        // save time of last render call
        vars->time_of_last_render_call_us = time_now_us;

        return PANE_HANDLER_RET_NO_ACTION;
    }

    // -----------------------
    // -------- EVENT --------
    // -----------------------

    if (request == PANE_HANDLER_REQ_EVENT) {
        if (IS_KEYBOARD_EVENT_ID(event->event_id)) {
            vars->key_event_id = event->event_id;
        } else {
            switch(event->event_id) {
            case SDL_EVENT_MOUSE_CLICK:
                vars->mouse_click_test_counter++;
                return PANE_HANDLER_RET_DISPLAY_REDRAW;
            case SDL_EVENT_MOUSE_WHEEL:
                vars->mouse_wheel_test_counter += event->mouse_wheel.delta_y;
                return PANE_HANDLER_RET_DISPLAY_REDRAW;
            case SDL_EVENT_NEW_DEMO_PANE:
                sdl_pane_create(pane_cx->pane_list_head, pane_handler_demo, NULL,
                                400, 400, pane_cx->w_total, pane_cx->h_total, PANE_BORDER_STYLE_STANDARD);
                return PANE_HANDLER_RET_DISPLAY_REDRAW;
            case SDL_EVENT_NEW_TEXT_PANE:
                sdl_pane_create(pane_cx->pane_list_head, pane_handler_display_text, text,
                                100, 100, 800, 800, PANE_BORDER_STYLE_STANDARD);
                return PANE_HANDLER_RET_DISPLAY_REDRAW;
            case SDL_EVENT_NEW_POINTS_PANE:
                sdl_pane_create(pane_cx->pane_list_head, pane_handler_points_test, NULL,
                                400, 0, 804, 804, PANE_BORDER_STYLE_MINIMAL);
                return PANE_HANDLER_RET_DISPLAY_REDRAW;
            case SDL_EVENT_NEW_DISPLAY:
                sdl_pane_manager(display_redraw_needed, 1,
                                 pane_handler_demo, NULL, 0, 0, 400, 400, PANE_BORDER_STYLE_STANDARD);
                return PANE_HANDLER_RET_DISPLAY_REDRAW;
            }
        }
        return PANE_HANDLER_RET_NO_ACTION;
    }

    // ---------------------------
    // -------- TERMINATE --------
    // ---------------------------

    if (request == PANE_HANDLER_REQ_TERMINATE) {
        sdl_destroy_texture(vars->circle_texture);
        free(vars);
        return PANE_HANDLER_RET_NO_ACTION;
    }

    // not reached
    assert(0);
    return PANE_HANDLER_RET_NO_ACTION;
}

static int32_t pane_handler_points_test(pane_cx_t * pane_cx, int32_t request, void * init, sdl_event_t * event) 
{
    #define SDL_EVENT_SELECT_POINT_SIZE     (SDL_EVENT_USER_DEFINED+0)

    #define MAX_POINTS 1000

    struct {
        point_t points[MAX_POINTS];
        int32_t point_size;
    } * vars = pane_cx->vars;

    // ----------------------------
    // -------- INITIALIZE --------
    // ----------------------------

    if (request == PANE_HANDLER_REQ_INITIALIZE) {
        int32_t i;
        vars = pane_cx->vars = calloc(1,sizeof(*vars));
        for (i = 0; i < MAX_POINTS; i++) {
            vars->points[i].x = random_range(0,pane_cx->pane.w-1);
            vars->points[i].y = random_range(0,pane_cx->pane.h-1);
        }
        vars->point_size = 0;
        return PANE_HANDLER_RET_NO_ACTION;
    }

    // ------------------------
    // -------- RENDER --------
    // ------------------------

    if (request == PANE_HANDLER_REQ_RENDER) {
        rect_t * pane = &pane_cx->pane;

        sdl_render_points(pane, vars->points, MAX_POINTS, BLUE, vars->point_size);

        sdl_render_text_and_register_event(
            pane, 0, 0, 0, "POINT_SIZE", LIGHT_BLUE, BLACK, 
            SDL_EVENT_SELECT_POINT_SIZE, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);
        sdl_render_printf(pane, 0, 11, 0, WHITE, BLACK, "%d", vars->point_size);

        return PANE_HANDLER_RET_NO_ACTION;
    }

    // -----------------------
    // -------- EVENT --------
    // -----------------------

    if (request == PANE_HANDLER_REQ_EVENT) {
        switch(event->event_id) {
        case SDL_EVENT_SELECT_POINT_SIZE:
            vars->point_size = (vars->point_size + 1) % 3;
            return PANE_HANDLER_RET_DISPLAY_REDRAW;
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

static bool display_redraw_needed(uint64_t time_render_us)
{
    // XXX is there a better way than this
    return microsec_timer() - time_render_us > 30000;
}

#if 0  // TEMPLATE
static int32_t pane_handler_xxx(pane_cx_t * pane_cx, int32_t request, sdl_event_t * event) 
{
    struct {
    } * vars = pane_cx->vars;

    // ----------------------------
    // -------- INITIALIZE --------
    // ----------------------------

    if (request == PANE_HANDLER_REQ_INITIALIZE) {
        vars = pane_cx->vars = calloc(1,sizeof(*vars));
        xxx
        return PANE_HANDLER_RET_NO_ACTION;
    }

    // ------------------------
    // -------- RENDER --------
    // ------------------------

    if (request == PANE_HANDLER_REQ_RENDER) {
        rect_t * pane = &pane_cx->pane;
        xxx
        return PANE_HANDLER_RET_NO_ACTION;
    }

    // -----------------------
    // -------- EVENT --------
    // -----------------------

    if (request == PANE_HANDLER_REQ_EVENT) {
        switch(event->event_id) {
        case SDL_EVENT_xxx
            return PANE_HANDLER_RET_DISPLAY_REDRAW;
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
#endif
