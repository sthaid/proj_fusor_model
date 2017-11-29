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

#include <pthread.h>
#include <math.h>

#include "model.h"
#include "util_sdl.h"
#include "util_misc.h"

//
// defines
//

#define DEFAULT_WIN_WIDTH  1900
#define DEFAULT_WIN_HEIGHT 1000

//
// typedefs
//

//
// variables
//

//
// prototypes
//

static void help(void);
static int32_t pane_hndlr_chamber(pane_cx_t * pane_cx, int32_t request, void * init, sdl_event_t * event) ;
static int32_t pane_hndlr_params(pane_cx_t * pane_cx, int32_t request, void * init, sdl_event_t * event) ;
static bool display_redraw_needed(uint64_t time_render_us);

// -----------------  MAIN  -------------------------------------------------------------

int main(int argc, char **argv)
{
    float chamber_radius   = 0.076;      // 3 inches
    float grid_radius      = 0.019;      // .75 inches
    float chamber_pressure = 3.9996711;  // 30 mTorr
    float grid_voltage     = -30000.0;   // - 30 kV
    float grid_current     = 0.005;      // 5 mA

    int32_t win_width, win_height; 

    // get and process options
    // -R <meters>  : chamber radius
    // -r <meters>  : grid radius
    // -p <pascals> : chamber pressure
    // -v <volts>   : grid voltage, must be negative
    // -c <amps>    : grid current
    // -h           : help
    while (true) {
        char opt_char = getopt(argc, argv, "R:r:p:v:c:h");
        if (opt_char == -1) {
            break;
        }
        switch (opt_char) {
        case 'R':
            if (sscanf(optarg, "%f", &chamber_radius) != 1) {
                FATAL("opt %c, '%s' not a number\n", opt_char, optarg);
            }
            break;
        case 'r':
            if (sscanf(optarg, "%f", &grid_radius) != 1) {
                FATAL("opt %c, '%s' not a number\n", opt_char, optarg);
            }
            break;
        case 'p':
            if (sscanf(optarg, "%f", &chamber_pressure) != 1) {
                FATAL("opt %c, '%s' not a number\n", opt_char, optarg);
            }
            break;
        case 'v':
            if (sscanf(optarg, "%f", &grid_voltage) != 1) {
                FATAL("opt %c, '%s' not a number\n", opt_char, optarg);
            }
            break;
        case 'c':
            if (sscanf(optarg, "%f", &grid_current) != 1) {
                FATAL("opt %c, '%s' not a number\n", opt_char, optarg);
            }
            break;
        case 'h':
            help();
            return 0;
        default:
            return 1;
            break;
        }
    }

    // initialize model
    model_init(chamber_radius, grid_radius, chamber_pressure, grid_voltage, grid_current);

    // use sdl to display the model
    win_width  = DEFAULT_WIN_WIDTH;
    win_height = DEFAULT_WIN_HEIGHT;
    if (sdl_init(&win_width, &win_height, true, false) < 0) {
        FATAL("sdl_init %dx%d failed\n", win_width, win_height);
    }
    sdl_pane_manager(
        display_redraw_needed, 2,
        pane_hndlr_chamber, NULL, 0,   0, 800, 800, PANE_BORDER_STYLE_MINIMAL,
        pane_hndlr_params,  NULL, 0, 800, 800, 200, PANE_BORDER_STYLE_MINIMAL);

    // terminate the model
    model_terminate();

    // terminate 
    return 0;
}

static void help(void)
{
    // xxx tbd
    printf("TBD\n");
}

// -----------------  PANE HANDLERS ----------------------------------------

static int32_t pane_hndlr_chamber(pane_cx_t * pane_cx, int32_t request, void * init, sdl_event_t * event) 
{
    #define SDL_EVENT_MOUSE_MOTION (SDL_EVENT_USER_DEFINED+0)
    #define SDL_EVENT_MOUSE_WHEEL  (SDL_EVENT_USER_DEFINED+1)
    #define SDL_EVENT_DISP_SELECT  (SDL_EVENT_USER_DEFINED+2)
    #define SDL_EVENT_STATE_SELECT (SDL_EVENT_USER_DEFINED+3)
    #define SDL_EVENT_GRID_SELECT  (SDL_EVENT_USER_DEFINED+4)

    enum disp {DISP_XY, DISP_XZ, DISP_YZ};

    #define STATE_STR \
        (model_is_running() ? "RUNNING" : "STOPPED")
    #define DISP_STR \
        (vars->disp == DISP_XY ? "XY" : vars->disp == DISP_XZ ? "XZ" : vars->disp == DISP_YZ ? "YZ"  : "??")

    struct {
        float x_offset;
        float y_offset;
        float width;
        float meters_per_pixel;
        int64_t updates;
        bool grid;
        enum disp disp;
        uint64_t time_us_last;
        double time_secs_last;
        float metric;
    } * vars = pane_cx->vars;
    rect_t * pane = &pane_cx->pane;

    // ----------------------------
    // -------- INITIALIZE --------
    // ----------------------------

    if (request == PANE_HANDLER_REQ_INITIALIZE) {
        vars = pane_cx->vars = calloc(1,sizeof(*vars));
        vars->x_offset = 0;
        vars->y_offset = 0;
        vars->width = params.chamber_radius * 2;
        vars->meters_per_pixel = vars->width / pane->w;
        vars->updates = 0;
        vars->grid = false;
        vars->disp = DISP_XY;
        vars->time_us_last = 0;
        vars->time_secs_last = 0;
        vars->metric = 0;
        return PANE_HANDLER_RET_NO_ACTION;
    }

    // ------------------------
    // -------- RENDER --------
    // ------------------------

    if (request == PANE_HANDLER_REQ_RENDER) {

        // render the points; use vars->disp to select the orientation
        {
        #define MAX_POINTS 1000
        #define POINT_SIZE 1
        #define ATOM_COLOR BLUE
        #define ION_COLOR  RED

        int32_t      idx1, idx2, x, y;
        particle_t * p;
        cell_t   *   c;
        point_t      points_atom[MAX_POINTS];
        point_t      points_ion[MAX_POINTS];
        int32_t      max_points_atom = 0, max_points_ion = 0;

        for (idx1 = 0; idx1 < MAX_CELL; idx1++) {
            for (idx2 = 0; idx2 < MAX_CELL; idx2++) {

                // xxx max_cell must be even  ??
                // xxx also c+1
                // xxx or don't bother with different orientations
                c = (vars->disp == DISP_XY ? &cell[idx1][idx2][MAX_CELL/2] :
                     vars->disp == DISP_XZ ? &cell[idx1][MAX_CELL/2][idx2] :
                                             &cell[MAX_CELL/2][idx1][idx2]);

                // xxx this needs locking
                LIST_FOREACH(p, &c->particle_list_head, cell_entries) {
                    if (vars->disp == DISP_XY) {
                        x = pane->w/2 + (p->x + vars->x_offset) / vars->meters_per_pixel; 
                        y = pane->h/2 + (p->y + vars->y_offset) / vars->meters_per_pixel;
                        // xxx assert(p->z >= 0 && p->z < CELL_SIZE);
                    } else if (vars->disp == DISP_XZ) {
                        x = pane->w/2 + (p->x + vars->x_offset) / vars->meters_per_pixel; 
                        y = pane->h/2 + (p->z + vars->y_offset) / vars->meters_per_pixel;
                        assert(p->y >= 0 && p->y < CELL_SIZE);
                    } else {  // vars->disp == YZ
                        x = pane->w/2 + (p->y + vars->x_offset) / vars->meters_per_pixel; 
                        y = pane->h/2 + (p->z + vars->y_offset) / vars->meters_per_pixel;
                        assert(p->x >= 0 && p->x < CELL_SIZE);
                    }

                    if (!p->ion) {
                        points_atom[max_points_atom].x = x;
                        points_atom[max_points_atom].y = y;
                        max_points_atom++;
                        if (max_points_atom == MAX_POINTS) {
                            sdl_render_points(pane, points_atom, max_points_atom, ATOM_COLOR, POINT_SIZE);
                            max_points_atom = 0;
                        }
                    } else {
                        points_ion[max_points_ion].x = x;
                        points_ion[max_points_ion].y = y;
                        max_points_ion++;
                        if (max_points_ion == MAX_POINTS) {
                            sdl_render_points(pane, points_ion, max_points_ion, ION_COLOR, POINT_SIZE);
                            max_points_ion = 0;
                        }
                    }
                }
            }
        }
        if (max_points_atom > 0) {
            sdl_render_points(pane, points_atom, max_points_atom, ATOM_COLOR, POINT_SIZE);
            max_points_atom = 0;
        }
        if (max_points_ion > 0) {
            sdl_render_points(pane, points_ion, max_points_ion, ION_COLOR, POINT_SIZE);
            max_points_ion = 0;
        }
        }

        // render the grid
        {
        float x, y, tmp, r, r_squared;
        r = params.chamber_radius;
        r_squared = r * r;
        if (vars->grid) {
            for (x = -r; x <= r; x += CELL_SIZE) {
                int32_t xpane, ypane_min, ypane_max;
                tmp = r_squared - x*x;
                if (tmp < 0) continue;
                tmp = sqrt(tmp);
                ypane_min = pane->h/2 + (vars->y_offset - tmp) / vars->meters_per_pixel;
                ypane_max = pane->h/2 + (vars->y_offset + tmp) / vars->meters_per_pixel;
                if (ypane_min < 0) ypane_min = 0;
                if (ypane_max > pane->h-1) ypane_max = pane->h-1;
                xpane = pane->w/2 + (x + vars->x_offset) / vars->meters_per_pixel; 
                sdl_render_line(pane, xpane, ypane_min, xpane, ypane_max, GRAY);
            }
            for (y = -r; y <= r; y += CELL_SIZE) {
                int32_t ypane, xpane_min, xpane_max;
                tmp = r_squared - y*y;
                if (tmp < 0) continue;
                tmp = sqrt(tmp);
                xpane_min = pane->w/2 + (vars->x_offset - tmp) / vars->meters_per_pixel;
                xpane_max = pane->w/2 + (vars->x_offset + tmp) / vars->meters_per_pixel;
                if (xpane_min < 0) xpane_min = 0;
                if (xpane_max > pane->w-1) xpane_max = pane->w-1;
                ypane = pane->h/2 + (y + vars->y_offset) / vars->meters_per_pixel; 
                sdl_render_line(pane, xpane_min, ypane, xpane_max, ypane, GRAY);
            }
        }
        }

        // render the chamber and grid
        {
        sdl_render_circle(pane, 
                          pane->w/2 + vars->x_offset / vars->meters_per_pixel,
                          pane->h/2 + vars->y_offset / vars->meters_per_pixel,
                          params.chamber_radius / vars->meters_per_pixel, 2, WHITE);

        sdl_render_circle(pane, 
                          pane->w/2 + vars->x_offset / vars->meters_per_pixel,
                          pane->h/2 + vars->y_offset / vars->meters_per_pixel,
                          params.grid_radius / vars->meters_per_pixel, 2, WHITE);
        }

        // xxx metric
        uint64_t time_us_now = microsec_timer();
        if (model_is_running()) {
            if (time_us_now - vars->time_us_last > 1000000) {
                if (vars->time_us_last != 0) {
                    vars->metric = ((time_secs - vars->time_secs_last) * 1e9) / 
                                   ((time_us_now - vars->time_us_last) / 1e6);
                }
                vars->time_us_last = time_us_now;
                vars->time_secs_last = time_secs;
            }
        } else {
            vars->time_us_last = time_us_now;
            vars->time_secs_last = time_secs;
            vars->metric = 0;
        }


        // render the controls and status
        {
        rect_t  locf = {0,0,pane->w,pane->h};
        sdl_register_event(pane, &locf, SDL_EVENT_MOUSE_MOTION, SDL_EVENT_TYPE_MOUSE_MOTION, pane_cx);
        sdl_register_event(pane, &locf, SDL_EVENT_MOUSE_WHEEL, SDL_EVENT_TYPE_MOUSE_WHEEL, pane_cx);

        sdl_render_printf(pane, 0, -9, 0, WHITE, BLACK, 
            "%9ld", ++vars->updates);
        sdl_render_printf(pane, 1, -9, 0, WHITE, BLACK, 
            "%9lf", vars->metric);

        sdl_render_printf(pane, 0, 0, 0, WHITE, BLACK, 
            "T = %.3lf us", time_secs*1000000.);
        sdl_render_printf(pane, 1, 0, 0, WHITE, BLACK, 
            "W = %0.3f m", vars->width); 
        sdl_render_text_and_register_event(pane, 2, 0, 0, STATE_STR, LIGHT_BLUE, BLACK, 
            SDL_EVENT_STATE_SELECT, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);
        sdl_render_text_and_register_event(pane, 3, 0, 0, "GRID", LIGHT_BLUE, BLACK, 
            SDL_EVENT_GRID_SELECT, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);
        sdl_render_text_and_register_event(pane, 4, 0, 0, DISP_STR, LIGHT_BLUE, BLACK, 
            SDL_EVENT_DISP_SELECT, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);
        }

        return PANE_HANDLER_RET_NO_ACTION;
    }

    // -----------------------
    // -------- EVENT --------
    // -----------------------

    if (request == PANE_HANDLER_REQ_EVENT) {
        switch(event->event_id) {
        case SDL_EVENT_MOUSE_MOTION: {
            vars->x_offset += event->mouse_motion.delta_x * vars->meters_per_pixel;
            vars->y_offset += event->mouse_motion.delta_y * vars->meters_per_pixel;
            return PANE_HANDLER_RET_DISPLAY_REDRAW; }
        case SDL_EVENT_MOUSE_WHEEL:
            if (event->mouse_wheel.delta_y > 0) {
                vars->width /= 1.1;
            } else if (event->mouse_wheel.delta_y < 0) {
                vars->width *= 1.1;
            }
            vars->meters_per_pixel = vars->width / pane->w;
            return PANE_HANDLER_RET_DISPLAY_REDRAW;
        case SDL_EVENT_DISP_SELECT:
            vars->disp = (vars->disp + 1) % 3;
            return PANE_HANDLER_RET_DISPLAY_REDRAW;
        case SDL_EVENT_STATE_SELECT:
            if (model_is_running()) {
                model_stop();
            } else {
                model_start();
            }
            return PANE_HANDLER_RET_DISPLAY_REDRAW;
        case SDL_EVENT_GRID_SELECT:
            vars->grid = !vars->grid;
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

static int32_t pane_hndlr_params(pane_cx_t * pane_cx, int32_t request, void * init, sdl_event_t * event) 
{
    struct {
        int32_t none_yet;
    } * vars = pane_cx->vars;
    rect_t * pane = &pane_cx->pane;

    // ----------------------------
    // -------- INITIALIZE --------
    // ----------------------------

    if (request == PANE_HANDLER_REQ_INITIALIZE) {
        vars = pane_cx->vars = calloc(1,sizeof(*vars));
        return PANE_HANDLER_RET_NO_ACTION;
    }

    // ------------------------
    // -------- RENDER --------
    // ------------------------

    if (request == PANE_HANDLER_REQ_RENDER) {
        sdl_render_printf(pane, 0, 0, 0, WHITE, BLACK, 
            "chamber_radius   = %0.3f m", params.chamber_radius);
        sdl_render_printf(pane, 1, 0, 0, WHITE, BLACK, 
            "grid_radius      = %0.3f m", params.grid_radius);
        sdl_render_printf(pane, 2, 0, 0, WHITE, BLACK, 
            "chamber_pressure = %0.3f pascal", params.chamber_pressure);
        sdl_render_printf(pane, 3, 0, 0, WHITE, BLACK, 
            "grid_voltage     = %0.0f V", params.grid_voltage);
        sdl_render_printf(pane, 4, 0, 0, WHITE, BLACK, 
            "grid_current     = %0.3f A", params.grid_current);
        return PANE_HANDLER_RET_NO_ACTION;
    }

    // -----------------------
    // -------- EVENT --------
    // -----------------------

    if (request == PANE_HANDLER_REQ_EVENT) {
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

// XXX pane_handler_graphs




static bool display_redraw_needed(uint64_t time_render_us)  // xxx or time_last_render_us
{
    uint64_t time_now = microsec_timer();

    return time_now - time_render_us > 100000;
}
