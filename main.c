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

#include "model.h"
#include "util_sdl.h"
#include "util_misc.h"

#if 0 // XXX TODO
HIGH
- add a control to adjust the pancake height, and center the pancake 
- model the ionization rate, and dipslay graph of the result
- display temperature graph, and test by setting the shell to high temperature
- when zoomed in then display the locbox grid
- xy xz yz
- display performance metric
- ways to increase performance, such as skipping locboxs that dont have ions, for a bit
- add electric force to locbox, and test by creating some ions,  say 1 percent
- use memory mapped file, OR try periodically writing the file
- add controls to start and stop the model
- use a scale table, with factor the eighth root of 2
MEDIUM
- change param values by left or right click OR mouse wheel
LOW
- try using SDL_RenderSetLogicalSize
#endif

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
static bool display_redraw_needed(uint64_t time_render_us);  // xxx make this optional

// -----------------  MAIN  -------------------------------------------------------------

int main(int argc, char **argv)
{
    static char default_params_str[] = "150,38,15,30,5";
    char * params_str = NULL;
    char * filename_str = NULL;
    int32_t win_width, win_height; 

    // get and process options
    // -p <params_str>   : example "150,38,15,30,5" 
    //                       150 = chamber diameter in mm
    //                       38  = grid diameter in mm
    //                       15  = pressure in mtorr
    //                       30  = voltage in kv
    //                       5   = current in ma
    // -f <filename_str> : saved model filename
    // -h                : help
    while (true) {
        static bool opt_p_supplied = false;
        static bool opt_f_supplied = false;
        char opt_char = getopt(argc, argv, "p:f:h");
        if (opt_char == -1) {
            break;
        }
        switch (opt_char) {
        case 'p':
            if (opt_p_supplied || opt_f_supplied) {
                FATAL("-p and -f can not be combined\n");
            }
            params_str = optarg;
            opt_p_supplied = true;
            break;
        case 'f':
            if (opt_p_supplied || opt_f_supplied) {
                FATAL("-p and -f can not be combined\n");
            }
            filename_str = optarg;
            opt_f_supplied = true;
            break;
        case 'h':
            help();
            exit(0);
        default:
            exit(1);
            break;
        }
    }

    // if neither params_str or filename supplied then
    // use default params_str
    if (!params_str && !filename_str) {
        params_str = default_params_str;
    }
 
    // initialize model
    if (params_str) {
        model_init_from_params(params_str);
    } else {
        model_init_from_file(filename_str);
    }

    // xxx 
    model_start();

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

    struct {
        int64_t x_offset_nm;
        int64_t y_offset_nm;
        int64_t width_nm;
        int64_t updates;  // xxx
    } * vars = pane_cx->vars;
    rect_t * pane = &pane_cx->pane;

    // ----------------------------
    // -------- INITIALIZE --------
    // ----------------------------

    if (request == PANE_HANDLER_REQ_INITIALIZE) {
        vars = pane_cx->vars = calloc(1,sizeof(*vars));
        vars->x_offset_nm = 0;
        vars->y_offset_nm = 0;
        vars->width_nm = params.chamber_radius_nm * 2;
        vars->updates = 0;
        return PANE_HANDLER_RET_NO_ACTION;
    }

    // ------------------------
    // -------- RENDER --------
    // ------------------------

    if (request == PANE_HANDLER_REQ_RENDER) {
        #define MAX_POINTS 1000
        #define POINT_SIZE 1
        #define ATOM_COLOR BLUE
        #define ION_COLOR  RED

        rect_t      locf = {0,0,pane->w,pane->h};
        int64_t      x_offset_nm = vars->x_offset_nm;
        int64_t      y_offset_nm = vars->y_offset_nm;
        int64_t      nm_per_pixel = vars->width_nm / pane->w;
        int32_t      x_idx, y_idx, x, y;
        particle_t * p;
        locbox_t   * lb;
        point_t      points_atom[MAX_POINTS];
        point_t      points_ion[MAX_POINTS];
        int32_t      max_points_atom = 0, max_points_ion = 0;

        for (x_idx = 0; x_idx < MAX_LOCBOX; x_idx++) {
            for (y_idx = 0; y_idx < MAX_LOCBOX; y_idx++) {
                lb = &locbox[x_idx][y_idx][MAX_LOCBOX/2];
                pthread_spin_lock(&lb->particle_list_spinlock);
                LIST_FOREACH(p, &lb->particle_list_head, entries) {
                    assert(p->z_nm >= 0 && p->z_nm < LOCBOX_SIZE_MM*1000000);
                    x = pane->w/2 + (p->x_nm + x_offset_nm) / nm_per_pixel; 
                    y = pane->h/2 + (p->y_nm + y_offset_nm) / nm_per_pixel;
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
                pthread_spin_unlock(&lb->particle_list_spinlock);
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

        sdl_render_circle(pane, 
                          pane->w/2 + x_offset_nm / nm_per_pixel,
                          pane->h/2 + y_offset_nm / nm_per_pixel,
                          params.chamber_radius_nm / nm_per_pixel, 2, WHITE);

        sdl_render_circle(pane, 
                          pane->w/2 + x_offset_nm / nm_per_pixel,
                          pane->h/2 + y_offset_nm / nm_per_pixel,
                          params.grid_radius_nm / nm_per_pixel, 2, WHITE);

        vars->updates++;
        sdl_render_printf(pane, 0, 0, 0, WHITE, BLACK, 
            "T = %.3lf us", time_ns/1000.);
        sdl_render_printf(pane, 1, 0, 0, WHITE, BLACK, 
            "W = %ld mm", vars->width_nm/1000000);
        sdl_render_printf(pane, 0, -9, 0, WHITE, BLACK, 
            "%9ld", vars->updates);

        sdl_register_event(pane, &locf, SDL_EVENT_MOUSE_MOTION, SDL_EVENT_TYPE_MOUSE_MOTION, pane_cx);
        sdl_register_event(pane, &locf, SDL_EVENT_MOUSE_WHEEL, SDL_EVENT_TYPE_MOUSE_WHEEL, pane_cx);

        return PANE_HANDLER_RET_NO_ACTION;
    }

    // -----------------------
    // -------- EVENT --------
    // -----------------------

    if (request == PANE_HANDLER_REQ_EVENT) {
        switch(event->event_id) {
        case SDL_EVENT_MOUSE_MOTION: {
            int64_t nm_per_pixel = vars->width_nm / pane->w;
            vars->x_offset_nm += event->mouse_motion.delta_x * nm_per_pixel;
            vars->y_offset_nm += event->mouse_motion.delta_y * nm_per_pixel;
            return PANE_HANDLER_RET_DISPLAY_REDRAW; }
        case SDL_EVENT_MOUSE_WHEEL:
            if (event->mouse_wheel.delta_y > 0) {
                vars->width_nm *= 1.1;
            } else if (event->mouse_wheel.delta_y < 0) {
                vars->width_nm /= 1.1;
            }
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
            "ChamberDiameter = %d mm", params.chamber_radius_nm*2/1000000);
        sdl_render_printf(pane, 1, 0, 0, WHITE, BLACK, 
            "GridDiameter    = %d mm", params.grid_radius_nm*2/1000000);
        sdl_render_printf(pane, 2, 0, 0, WHITE, BLACK, 
            "ChamberPressure = %.1lf mTorr", params.chamber_pressure_utorr/1000.);
        sdl_render_printf(pane, 3, 0, 0, WHITE, BLACK, 
            "GridVoltage     = %.1lf kV", params.grid_voltage_v/1000.);
        sdl_render_printf(pane, 4, 0, 0, WHITE, BLACK, 
            "GridCurrent     = %.1lf mA", params.grid_current_ua/1000.);
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

static bool display_redraw_needed(uint64_t time_render_us)
{
    uint64_t time_now = microsec_timer();

    return time_now - time_render_us > 1000000;
}
