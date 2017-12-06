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

#include "util_sdl.h"
#include "util_sdl_predefined_panes.h"
#include "util_misc.h"

#include "model.h"
#include "physics.h"

//
// defines
//

#define DEFAULT_WIN_WIDTH  1900
#define DEFAULT_WIN_HEIGHT 1000

//
// typedefs
//

typedef struct {
    int32_t max;
    int32_t max_alloced;
    struct loc_s {
        float x;
        float y;
    } loc[0];
} loc_list_t;

//
// variables
//

loc_list_t                        * loc_atoms;
loc_list_t                        * loc_ions;
pane_hndlr_display_graph_params_t * gr_temperature;
pane_hndlr_display_graph_params_t * gr_nd_atoms;    

//
// prototypes
//

static void help(void);
static int32_t pane_hndlr_chamber(pane_cx_t * pane_cx, int32_t request, void * init_params, sdl_event_t * event) ;
static int32_t pane_hndlr_params(pane_cx_t * pane_cx, int32_t request, void * init_params, sdl_event_t * event) ;
static void display_start(void * display_cx);
static void display_end(void * display_cx);

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
        NULL,           // context
        display_start,  // called prior to pane handlers
        display_end,    // called after pane handlers
        250000,         // 0=continuous, -1=never, else us 
        4,              // number of pane handler varargs that follow
        pane_hndlr_chamber, NULL,     0,   0, 800, 800, PANE_BORDER_STYLE_MINIMAL,
        pane_hndlr_params,  NULL,     0, 800, 800, 200, PANE_BORDER_STYLE_MINIMAL,
        pane_hndlr_display_graph,  &gr_temperature, 800,   0, 550, 200, PANE_BORDER_STYLE_MINIMAL,
        pane_hndlr_display_graph,  &gr_nd_atoms,    800, 200, 550, 200, PANE_BORDER_STYLE_MINIMAL);

    // terminate 
    return 0;
}

static void help(void)
{
    // xxx tbd
    printf("TBD\n");
}

// -----------------  PANE HANDLERS ----------------------------------------

static int32_t pane_hndlr_chamber(pane_cx_t * pane_cx, int32_t request, void * init_params, sdl_event_t * event) 
{
    #define SDL_EVENT_MOUSE_MOTION (SDL_EVENT_USER_DEFINED+0)
    #define SDL_EVENT_MOUSE_WHEEL  (SDL_EVENT_USER_DEFINED+1)
    #define SDL_EVENT_STATE_SELECT (SDL_EVENT_USER_DEFINED+2)
    #define SDL_EVENT_STEP_SELECT  (SDL_EVENT_USER_DEFINED+3)
    #define SDL_EVENT_GRID_SELECT  (SDL_EVENT_USER_DEFINED+4)
    #define SDL_EVENT_RESET        (SDL_EVENT_USER_DEFINED+5)

    #define STATE_STR \
        (model_is_running() ? "RUNNING" : "STOPPED")

    struct {
        float x_offset;
        float y_offset;
        float width;
        float meters_per_pixel;
        bool grid;
        double time_wall_secs_last;
        double time_model_secs_last;
        float model_secs_per_wall_sec;
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
        vars->grid = false;
        vars->time_wall_secs_last = 0;
        vars->time_model_secs_last = 0;
        vars->model_secs_per_wall_sec = 0;
        return PANE_HANDLER_RET_NO_ACTION;
    }

    // ------------------------
    // -------- RENDER --------
    // ------------------------

    if (request == PANE_HANDLER_REQ_RENDER) {

        // render the location of the atoms & ions in the chamber 
        // which are located near z equals 0
        {
        #define MAX_POINTS 1000
        #define POINT_SIZE 1
        #define ATOM_COLOR BLUE
        #define ION_COLOR  RED

        point_t points[MAX_POINTS];
        int32_t i, max_points;

        max_points = 0;
        for (i = 0; i < loc_atoms->max; i++) {
            points[max_points].x = pane->w/2 + (loc_atoms->loc[i].x + vars->x_offset) / vars->meters_per_pixel; 
            points[max_points].y = pane->h/2 + (loc_atoms->loc[i].y + vars->y_offset) / vars->meters_per_pixel;
            max_points++;
            if (max_points == MAX_POINTS || i == loc_atoms->max-1) {
                sdl_render_points(pane, points, max_points, ATOM_COLOR, POINT_SIZE);
                max_points = 0;
            }
        }
        max_points = 0;
        for (i = 0; i < loc_ions->max; i++) {
            points[max_points].x = pane->w/2 + (loc_ions->loc[i].x + vars->x_offset) / vars->meters_per_pixel; 
            points[max_points].y = pane->h/2 + (loc_ions->loc[i].y + vars->y_offset) / vars->meters_per_pixel;
            max_points++;
            if (max_points == MAX_POINTS || i == loc_ions->max-1) {
                sdl_render_points(pane, points, max_points, ATOM_COLOR, POINT_SIZE);
                max_points = 0;
            }
        }
        }

        // if the grid enabled hen render the grid
        {
        float x, y, tmp, r, r_squared;
        r = params.chamber_radius;
        r_squared = r * r;
        if (vars->grid) {
            for (x = -r; x <= r; x += CELL_SIZE) {
                int32_t xpane, ypane_min, ypane_max;
                tmp = r_squared - x*x;
                if (tmp < 0) continue;
                tmp = sqrtf(tmp);
                ypane_min = pane->h/2 + (vars->y_offset - tmp) / vars->meters_per_pixel;
                ypane_max = pane->h/2 + (vars->y_offset + tmp) / vars->meters_per_pixel;
                if (ypane_min < 0) ypane_min = 0;
                if (ypane_max > pane->h-1) ypane_max = pane->h-1;
                xpane = pane->w/2 + (x + vars->x_offset) / vars->meters_per_pixel; 
                sdl_render_line(pane, xpane, ypane_min, xpane, ypane_max, GRAY);
                if (fabs(x) < CELL_SIZE/10) {
                    ypane_min = pane->h/2 + (vars->y_offset - CELL_SIZE) / vars->meters_per_pixel;
                    ypane_max = pane->h/2 + (vars->y_offset + CELL_SIZE) / vars->meters_per_pixel;
                    if (ypane_min < 0) ypane_min = 0;
                    if (ypane_max > pane->h-1) ypane_max = pane->h-1;
                    sdl_render_line(pane, xpane, ypane_min, xpane, ypane_max, RED);
                }
            }
            for (y = -r; y <= r; y += CELL_SIZE) {
                int32_t ypane, xpane_min, xpane_max;
                tmp = r_squared - y*y;
                if (tmp < 0) continue;
                tmp = sqrtf(tmp);
                xpane_min = pane->w/2 + (vars->x_offset - tmp) / vars->meters_per_pixel;
                xpane_max = pane->w/2 + (vars->x_offset + tmp) / vars->meters_per_pixel;
                if (xpane_min < 0) xpane_min = 0;
                if (xpane_max > pane->w-1) xpane_max = pane->w-1;
                ypane = pane->h/2 + (y + vars->y_offset) / vars->meters_per_pixel; 
                sdl_render_line(pane, xpane_min, ypane, xpane_max, ypane, GRAY);
                if (fabs(y) < CELL_SIZE/10) {
                    xpane_min = pane->w/2 + (vars->x_offset - CELL_SIZE) / vars->meters_per_pixel;
                    xpane_max = pane->w/2 + (vars->x_offset + CELL_SIZE) / vars->meters_per_pixel;
                    if (xpane_min < 0) xpane_min = 0;
                    if (xpane_max > pane->w-1) xpane_max = pane->w-1;
                    sdl_render_line(pane, xpane_min, ypane, xpane_max, ypane, RED);
                }
            }
        }
        }

        // render circles to represent the chamber and grid
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

        // determine the model progress rate, model_secs_per_wall_sec
        if (model_is_running()) {
            double time_wall_secs = microsec_timer() / 1000000.;
            if (vars->time_wall_secs_last == 0) {
                vars->time_wall_secs_last = time_wall_secs;
                vars->time_model_secs_last = time_model_secs;
            } else if (time_wall_secs - vars->time_wall_secs_last > 1) {
                vars->model_secs_per_wall_sec = 
                    (time_model_secs - vars->time_model_secs_last)  / 
                    (time_wall_secs - vars->time_wall_secs_last);
                vars->time_wall_secs_last = time_wall_secs;
                vars->time_model_secs_last = time_model_secs;
            }
        } else {
            vars->time_wall_secs_last = 0;
            vars->time_model_secs_last = 0;
            vars->model_secs_per_wall_sec = 0;
        }

        // render the controls and status
        {
        char   str[100];
        rect_t locf = {0,0,pane->w,pane->h};
        // status
        sdl_render_printf(pane, COL2X(0,1), ROW2Y(0,1), 1, WHITE, BLACK, 
            "T = %.3f us", time_model_secs*1e6);
        sdl_render_printf(pane, COL2X(0,1), ROW2Y(1,1), 1, WHITE, BLACK, 
            "W = %.3f m", vars->width); 
        if (vars->model_secs_per_wall_sec != 0) {
            sprintf(str, "%.3f us/s", vars->model_secs_per_wall_sec*1e6);
            sdl_render_text(pane, COL2X(-strlen(str),1), ROW2Y(0,1), 1, str, WHITE, BLACK);
        }
        // controls
        sdl_register_event(pane, &locf, SDL_EVENT_MOUSE_MOTION, SDL_EVENT_TYPE_MOUSE_MOTION, pane_cx);
        sdl_register_event(pane, &locf, SDL_EVENT_MOUSE_WHEEL, SDL_EVENT_TYPE_MOUSE_WHEEL, pane_cx);
        sdl_render_text_and_register_event(pane, COL2X(0,1), ROW2Y(2,1), 1, STATE_STR, LIGHT_BLUE, BLACK, 
            SDL_EVENT_STATE_SELECT, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);
        if (!model_is_running() && !model_is_stepping()) {
            sdl_render_text_and_register_event(pane, COL2X(0,1), ROW2Y(3,1), 1, "STEP", LIGHT_BLUE, BLACK, 
                SDL_EVENT_STEP_SELECT, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);
        }
        sdl_render_text_and_register_event(pane, COL2X(0,1), ROW2Y(4,1), 1, "GRID", LIGHT_BLUE, BLACK, 
            SDL_EVENT_GRID_SELECT, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);
        sdl_render_text_and_register_event(pane, COL2X(0,1), ROW2Y(5,1), 1, "R", LIGHT_BLUE, BLACK, 
            SDL_EVENT_RESET, SDL_EVENT_TYPE_MOUSE_CLICK, pane_cx);
        }

        // return
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
        case SDL_EVENT_STATE_SELECT:
            if (model_is_running()) {
                model_stop();
            } else {
                model_run();
            }
            return PANE_HANDLER_RET_DISPLAY_REDRAW;
        case SDL_EVENT_STEP_SELECT:
            model_step();
            return PANE_HANDLER_RET_DISPLAY_REDRAW;
        case SDL_EVENT_GRID_SELECT:
            vars->grid = !vars->grid;
            return PANE_HANDLER_RET_DISPLAY_REDRAW;
        case SDL_EVENT_RESET:
            vars->x_offset = 0;
            vars->y_offset = 0;
            vars->width = params.chamber_radius * 2;
            vars->meters_per_pixel = vars->width / pane->w;
            vars->grid = false;
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

static int32_t pane_hndlr_params(pane_cx_t * pane_cx, int32_t request, void * init_params, sdl_event_t * event) 
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
        sdl_render_printf(pane, COL2X(0,1), ROW2Y(0,1), 1, WHITE, BLACK, 
            "chamber_radius   = %0.3f m", params.chamber_radius);
        sdl_render_printf(pane, COL2X(0,1), ROW2Y(1,1), 1, WHITE, BLACK, 
            "grid_radius      = %0.3f m", params.grid_radius);
        sdl_render_printf(pane, COL2X(0,1), ROW2Y(2,1), 1, WHITE, BLACK, 
            "chamber_pressure = %0.3f pascal", params.chamber_pressure);
        sdl_render_printf(pane, COL2X(0,1), ROW2Y(3,1), 1, WHITE, BLACK, 
            "grid_voltage     = %0.0f V", params.grid_voltage);
        sdl_render_printf(pane, COL2X(0,1), ROW2Y(4,1), 1, WHITE, BLACK, 
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

static void display_start(void * display_cx)
{
    // XXX overview comment

    // XXX check time since last< don't exceed

    // pause the model
    model_pause();

    // create list of atom/ion locations with z approx 0
    // XXX max_cell must be even  ??
    {
    int32_t x_idx, y_idx;
    cell_t * c;
    particle_t * p;
    if (loc_atoms == NULL) {
        loc_atoms = malloc(sizeof(*loc_atoms) + 1000000*sizeof(struct loc_s));
        loc_atoms->max = 0;
        loc_atoms->max_alloced = 1000000;
    }
    loc_atoms->max = 0;
    if (loc_ions == NULL) {
        loc_ions = malloc(sizeof(*loc_ions) + 1000000*sizeof(struct loc_s));
        loc_ions->max = 0;
        loc_ions->max_alloced = 1000000;
    }
    loc_ions->max = 0;
    for (x_idx = 0; x_idx < MAX_CELL; x_idx++) {
        for (y_idx = 0; y_idx < MAX_CELL; y_idx++) {
            c = &cell[x_idx][y_idx][MAX_CELL/2];
            LIST_FOREACH(p, &c->particle_list_head, cell_entries) {
                assert(p->z >= 0 && p->z < CELL_SIZE);
                if (!p->ion) {
                    loc_atoms->loc[loc_atoms->max].x = p->x;
                    loc_atoms->loc[loc_atoms->max].y = p->y;
                    loc_atoms->max++;
                    if (loc_atoms->max == loc_atoms->max_alloced) {
                        loc_atoms->max_alloced += 1000000;
                        loc_atoms = realloc(loc_atoms, 
                                            sizeof(*loc_atoms) + loc_atoms->max_alloced*sizeof(struct loc_s));
                        DEBUG("XXX realloced loc_atoms  %d\n", loc_atoms->max_alloced);
                    }
                } else {
                    loc_ions->loc[loc_ions->max].x = p->x;
                    loc_ions->loc[loc_ions->max].y = p->y;
                    loc_ions->max++;
                    if (loc_ions->max == loc_ions->max_alloced) {
                        loc_ions->max_alloced += 1000000;
                        loc_ions = realloc(loc_ions, 
                                            sizeof(*loc_ions) + loc_ions->max_alloced*sizeof(struct loc_s));
                        DEBUG("XXX realloced loc_ions  %d\n", loc_ions->max_alloced);
                    }
                }
            }
        }
    }
    }

    // graph x=radius y=temperature  
    {
    int32_t i, max_points_needed;
    bool init = false;
    max_points_needed = max_shell;
    if (gr_temperature == NULL || (gr_temperature)->max_points_alloced < max_points_needed) {
        free(gr_temperature);
        gr_temperature = 
            malloc(sizeof(pane_hndlr_display_graph_params_t) +
                   max_points_needed * sizeof(struct pane_hndlr_display_graph_point_s));
        init = true;
    }
    pane_hndlr_display_graph_params_t * g = gr_temperature;
    if (init) {
        strcpy(g->title_str, "TEMPERATURE");
        strcpy(g->x_units_str, "METERS");
        strcpy(g->y_units_str, "K");
        g->x_min = 0;
        g->x_max = max_shell * SHELL_SIZE;
        g->y_min = 0;
        g->y_max = 600;
        g->max_points_alloced = max_points_needed;
    }
    for (i = 0; i < max_points_needed; i++) {
        double avg_v_squared = shell[i].sum_v_squared / (shell[i].number_of_atoms + shell[i].number_of_ions);
        double temperature = V_SQUARED_TO_TEMPERATURE(avg_v_squared,D_MASS);
        g->points[i].x = i * SHELL_SIZE;
        g->points[i].y = temperature;
    }
    g->max_points = max_points_needed;
    }

    // graph x=radius y=nd_atoms  
    {
    int32_t i, max_points_needed;
    bool init = false;
    max_points_needed = max_shell;
    if (gr_nd_atoms == NULL || (gr_nd_atoms)->max_points_alloced < max_points_needed) {
        free(gr_nd_atoms);
        gr_nd_atoms = 
            malloc(sizeof(pane_hndlr_display_graph_params_t) +
                   max_points_needed * sizeof(struct pane_hndlr_display_graph_point_s));
        init = true;
    }
    pane_hndlr_display_graph_params_t * g = gr_nd_atoms;
    if (init) {
        strcpy(g->title_str, "ND_ATOMS");
        strcpy(g->x_units_str, "METERS");
        strcpy(g->y_units_str, "/M^3");
        g->x_min = 0;
        g->x_max = max_shell * SHELL_SIZE;
        g->y_min = 0;
        g->y_max = 2e21;  // xxx calc from pressure, etc
        g->max_points_alloced = max_points_needed;
    }
    for (i = 0; i < max_points_needed; i++) {
        g->points[i].x = i * SHELL_SIZE;
        g->points[i].y = shell[i].number_of_atoms * num_real_particles_per_sim_particle / shell[i].volume;
    }
    g->max_points = max_points_needed;
    }

    // resume the model
    model_resume();
}

static void display_end(void * display_cx)
{
}
