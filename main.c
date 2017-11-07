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

#if 0 // XXX TODO
HIGH
- display temperature graph, and test by setting the shell to high temperature
- ways to increase performance, such as skipping locboxs that dont have ions, for a bit
- add a control to adjust the pancake height, and center the pancake 
- model the ionization rate, and dipslay graph of the result
- display performance metric
- use memory mapped file, OR try periodically writing the file
MEDIUM
- change param values by left or right click OR mouse wheel
- use a scale table, with factor the eighth root of 2
LOW
- try using SDL_RenderSetLogicalSize

DONE
- xy xz yz
- add electric force to locbox, and test by creating some ions,  say 1 percent
- add controls to start and stop the model
- when zoomed in then display the locbox grid
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
    #define SDL_EVENT_DISP_SELECT  (SDL_EVENT_USER_DEFINED+2)
    #define SDL_EVENT_STATE_SELECT (SDL_EVENT_USER_DEFINED+3)
    #define SDL_EVENT_GRID_SELECT  (SDL_EVENT_USER_DEFINED+4)

    #define STATE_RUNNING 0
    #define STATE_STOPPED 1
    #define STATE_STR \
        (vars->state == STATE_RUNNING ? "RUNNING" : vars->state == STATE_STOPPED ? "STOPPED" : "??")

    #define GRID_OFF 0
    #define GRID_ON  1
    #define GRID_STR \
        "GRID"  // xxx GRID_ON or GRID_OFF  ????

    #define DISP_XY 0
    #define DISP_XZ 1
    #define DISP_YZ 2
    #define DISP_STR \
        (vars->disp == DISP_XY ? "XY" : vars->disp == DISP_XZ ? "XZ" : vars->disp == DISP_YZ ? "YZ"  : "??")

    #define NM_PER_PIXEL (vars->width_nm / pane->w)

    struct {
        int64_t x_offset_nm;
        int64_t y_offset_nm;
        int64_t width_nm;
        int64_t updates;  // xxx
        int64_t state;  // xxx bool or enum
        int64_t grid;  // xxx bool
        int64_t disp;   // xxx enum
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
        vars->state = STATE_RUNNING;   // xxx maybe should query the model
        vars->grid = GRID_OFF;
        vars->disp = DISP_XY;
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
        locbox_t   * lb;
        point_t      points_atom[MAX_POINTS];
        point_t      points_ion[MAX_POINTS];
        int32_t      max_points_atom = 0, max_points_ion = 0;

        for (idx1 = 0; idx1 < MAX_LOCBOX; idx1++) {
            for (idx2 = 0; idx2 < MAX_LOCBOX; idx2++) {

                lb = (vars->disp == DISP_XY ? &locbox[idx1][idx2][MAX_LOCBOX/2] :
                      vars->disp == DISP_XZ ? &locbox[idx1][MAX_LOCBOX/2][idx2] :
                                              &locbox[MAX_LOCBOX/2][idx1][idx2]);

                pthread_spin_lock(&lb->particle_list_spinlock);
                LIST_FOREACH(p, &lb->particle_list_head, entries) {
                    if (vars->disp == DISP_XY) {
                        x = pane->w/2 + (p->x_nm + vars->x_offset_nm) / NM_PER_PIXEL; 
                        y = pane->h/2 + (p->y_nm + vars->y_offset_nm) / NM_PER_PIXEL;
                        assert(p->z_nm >= 0 && p->z_nm < LOCBOX_SIZE_NM);
                    } else if (vars->disp == DISP_XZ) {
                        x = pane->w/2 + (p->x_nm + vars->x_offset_nm) / NM_PER_PIXEL; 
                        y = pane->h/2 + (p->z_nm + vars->y_offset_nm) / NM_PER_PIXEL;
                        assert(p->y_nm >= 0 && p->y_nm < LOCBOX_SIZE_NM);
                    } else {  // vars->disp == YZ
                        x = pane->w/2 + (p->y_nm + vars->x_offset_nm) / NM_PER_PIXEL; 
                        y = pane->h/2 + (p->z_nm + vars->y_offset_nm) / NM_PER_PIXEL;
                        assert(p->x_nm >= 0 && p->x_nm < LOCBOX_SIZE_NM);
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
        }

        // render the grid
        {
        int64_t x_nm, y_nm, tmp_nm, r_nm = params.chamber_radius_nm;
        int64_t x, y_min, y_max;
        int64_t y, x_min, x_max;
        if (vars->grid == GRID_ON) {
            for (x_nm = -r_nm; x_nm <= r_nm; x_nm += LOCBOX_SIZE_NM) {
                x = pane->w/2 + (x_nm + vars->x_offset_nm) / NM_PER_PIXEL; 
                tmp_nm = r_nm*r_nm - x_nm*x_nm;
                if (tmp_nm < 0) continue;
                tmp_nm = sqrt(tmp_nm);
                y_min = pane->h/2 + (vars->y_offset_nm - tmp_nm) / NM_PER_PIXEL;
                y_max = pane->h/2 + (vars->y_offset_nm + tmp_nm) / NM_PER_PIXEL;
                if (y_min < 0) y_min = 0;
                if (y_max > pane->h-1) y_max = pane->h-1;
                sdl_render_line(pane, x, y_min, x, y_max, GRAY);
            }
            for (y_nm = -r_nm; y_nm <= r_nm; y_nm += LOCBOX_SIZE_NM) {
                y = pane->h/2 + (y_nm + vars->y_offset_nm) / NM_PER_PIXEL; 
                tmp_nm = r_nm*r_nm - y_nm*y_nm;
                if (tmp_nm < 0) continue;
                tmp_nm = sqrt(tmp_nm);
                x_min = pane->w/2 + (vars->x_offset_nm - tmp_nm) / NM_PER_PIXEL;
                x_max = pane->w/2 + (vars->x_offset_nm + tmp_nm) / NM_PER_PIXEL;
                if (x_min < 0) x_min = 0;
                if (x_max > pane->w-1) x_max = pane->w-1;
                sdl_render_line(pane, x_min, y, x_max, y, GRAY);
            }
        }
        }

        // render the chamber and grid
        {
        sdl_render_circle(pane, 
                          pane->w/2 + vars->x_offset_nm / NM_PER_PIXEL,
                          pane->h/2 + vars->y_offset_nm / NM_PER_PIXEL,
                          params.chamber_radius_nm / NM_PER_PIXEL, 2, WHITE);

        sdl_render_circle(pane, 
                          pane->w/2 + vars->x_offset_nm / NM_PER_PIXEL,
                          pane->h/2 + vars->y_offset_nm / NM_PER_PIXEL,
                          params.grid_radius_nm / NM_PER_PIXEL, 2, WHITE);
        }

        // render the controls and status
        {
        rect_t  locf = {0,0,pane->w,pane->h};
        sdl_register_event(pane, &locf, SDL_EVENT_MOUSE_MOTION, SDL_EVENT_TYPE_MOUSE_MOTION, pane_cx);
        sdl_register_event(pane, &locf, SDL_EVENT_MOUSE_WHEEL, SDL_EVENT_TYPE_MOUSE_WHEEL, pane_cx);

        sdl_render_printf(pane, 0, -9, 0, WHITE, BLACK, 
            "%9ld", ++vars->updates);

        sdl_render_printf(pane, 0, 0, 0, WHITE, BLACK, 
            "T = %.3lf us", time_ns/1000.);
        sdl_render_printf(pane, 1, 0, 0, WHITE, BLACK, 
            "W = %ld mm", vars->width_nm/1000000);  // xxx conversion macro
        sdl_render_text_and_register_event(pane, 2, 0, 0, STATE_STR, LIGHT_BLUE, BLACK, 
            SDL_EVENT_STATE_SELECT, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);
        sdl_render_text_and_register_event(pane, 3, 0, 0, GRID_STR, LIGHT_BLUE, BLACK, 
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
            vars->x_offset_nm += event->mouse_motion.delta_x * NM_PER_PIXEL;
            vars->y_offset_nm += event->mouse_motion.delta_y * NM_PER_PIXEL;
            return PANE_HANDLER_RET_DISPLAY_REDRAW; }
        case SDL_EVENT_MOUSE_WHEEL:
            if (event->mouse_wheel.delta_y > 0) {
                vars->width_nm /= 1.1;
            } else if (event->mouse_wheel.delta_y < 0) {
                vars->width_nm *= 1.1;
            }
            return PANE_HANDLER_RET_DISPLAY_REDRAW;
        case SDL_EVENT_DISP_SELECT:
            vars->disp = (vars->disp + 1) % 3;
            return PANE_HANDLER_RET_DISPLAY_REDRAW;
        case SDL_EVENT_STATE_SELECT:
            if (vars->state == STATE_RUNNING) {
                vars->state = STATE_STOPPED;
                model_stop();
            } else {   // vars->state == STATE_STOPPED
                vars->state = STATE_RUNNING;
                model_start();
            }
            return PANE_HANDLER_RET_DISPLAY_REDRAW;
        case SDL_EVENT_GRID_SELECT:
            vars->grid = (vars->grid == GRID_ON ? GRID_OFF : GRID_ON);
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
            "ChamberDiameter = %d mm", params.chamber_radius_nm*2/1000000);  // xxx conversion macro
        sdl_render_printf(pane, 1, 0, 0, WHITE, BLACK, 
            "GridDiameter    = %d mm", params.grid_radius_nm*2/1000000);  // xxx conversion macro
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
