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

#include <math.h>
#include <float.h>
#include <pthread.h>
#include <sys/mman.h>

#include "util_sdl.h"
#include "util_sdl_predefined_panes.h"
#include "util_misc.h"

#include "physics.h"
#include "model.h"

//
// defines
//

// model time increment
#define MODEL_TIME_INCREMENT 1e-9   // 1 ns

// misc
#define MB 0x100000

// 
// typedefs
//

//
// variables
//

static bool run_request;
static volatile bool running;
static bool model_thread_started;
static float roomtemp_d_velocity;
static float roomtemp_d_velocity_squared;
 
//
// prototypes
//

static void init_particle(particle_t * p);
static void * model_thread(void * cx);

//
// inline functions
//

// xxx comment
static inline cell_t * get_cell(float x, float y, float z)
{
    cell_t * c;

    // compute the x/y/z_idx for the cell
    int32_t x_idx = (x + (MAX_CELL*CELL_SIZE/2)) / CELL_SIZE;
    int32_t y_idx = (y + (MAX_CELL*CELL_SIZE/2)) / CELL_SIZE;
    int32_t z_idx = (z + (MAX_CELL*CELL_SIZE/2)) / CELL_SIZE;

    // assert the idx values are in range
    assert(x_idx >= 0 && x_idx < MAX_CELL);
    assert(y_idx >= 0 && y_idx < MAX_CELL);
    assert(z_idx >= 0 && z_idx < MAX_CELL);

    // if cell is not in chamber then return null
    c = &cell[x_idx][y_idx][z_idx];
    if (c->shell == NULL) {
        return NULL;
    }
    
    // return ptr to cell
    return c;
}

// xxx comment
static inline void random_location_within_chamber(float *x, float *y, float *z, cell_t **c)
{
    float x_try, y_try, z_try;
    cell_t * c_try;

    while (true) {
        x_try = random_range(-params.chamber_radius, params.chamber_radius);
        y_try = random_range(-params.chamber_radius, params.chamber_radius);
        z_try = random_range(-params.chamber_radius, params.chamber_radius);
        c_try = get_cell(x_try, y_try, z_try);
        if (c_try) {
            *x = x_try;
            *y = y_try;
            *z = z_try;
            *c = c_try;
            break;
        }
    }
}

// -----------------  MODEL_INIT  --------------------------------------------------------

void model_init(float chamber_radius, float grid_radius, float chamber_pressure,
        float grid_voltage, float grid_current)
{
    int32_t ret;

    // debug print sizes
    DEBUG("MAX_PARTICLES     = %d\n", MAX_PARTICLES);
    DEBUG("MAX_CELL          = %d\n", MAX_CELL);
    DEBUG("MAX_SHELL         = %d\n", MAX_SHELL);
    DEBUG("sizeof(particles) = %ld MB\n", sizeof(particles) / MB);
    DEBUG("sizeof(cell)      = %ld MB\n", sizeof(cell) / MB);
    DEBUG("sizeof(shell)     = %ld MB\n", sizeof(shell) / MB);
    DEBUG("sizeof(double)    = %ld bytes\n", sizeof(double));
    DEBUG("sizeof(float)     = %ld bytes\n", sizeof(float));
    DEBUG("sizeof(1.23)      = %ld bytes\n", sizeof(1.23));
    DEBUG("DBL_MAX           = %e\n", DBL_MAX);
    DEBUG("DBL_MIN           = %e\n", DBL_MIN);
    DEBUG("FLT_MAX           = %e\n", FLT_MAX);
    DEBUG("FLT_MIN           = %e\n", FLT_MIN);

    // xxx comment and mcmodel
    // -mcmodel=medium
    particles = calloc(MAX_PARTICLES, sizeof(particle_t));
    if (particles == NULL) {
        FATAL("failed calloc particles, %ld MB\n", MAX_PARTICLES*sizeof(particle_t));
    }

    // don't dump the particles and cell arrays, because they are so large
    ret = madvise((void*)ROUND_UP(particles,PAGE_SIZE), 
                  MAX_PARTICLES*sizeof(particle_t)-PAGE_SIZE, 
                  MADV_DONTDUMP);
    if (ret != 0) {
        FATAL("madvise particles ret %d, %s\n", ret, strerror(errno));
    }
    ret = madvise((void*)ROUND_UP(cell,PAGE_SIZE), 
                  sizeof(cell)-PAGE_SIZE, 
                  MADV_DONTDUMP);
    if (ret != 0) {
        FATAL("madvise cell ret %d, %s\n", ret, strerror(errno));
    }

    // verify params are in range
    if (chamber_radius < MIN_CHAMBER_RADIUS ||
        chamber_radius > MAX_CHAMBER_RADIUS) {
        FATAL("chamber_radius %lf is out of range\n", chamber_radius);
    }
    if (grid_radius < MIN_GRID_RADIUS ||
        grid_radius > MAX_GRID_RADIUS(chamber_radius)) {
        FATAL("grid_radius %lf is out of range\n", grid_radius);
    }
    if (chamber_pressure < MIN_CHAMBER_PRESSURE ||
        chamber_pressure > MAX_CHAMBER_PRESSURE) {
        FATAL("chamber_pressure %lf is out of range\n", chamber_pressure);
    }
    if (grid_voltage > 0 || 
        grid_voltage < MAX_GRID_VOLTAGE)  {
        FATAL("grid_voltage %lf is out of range\n", grid_voltage);
    }
    if (grid_current < 0 ||
        grid_current > MAX_GRID_CURRENT) {
        FATAL("grid_current %lf is out of range\n", grid_current);
    }

    // copy the validated param values in to the params struct
    params.chamber_radius    = chamber_radius;
    params.grid_radius       = grid_radius;
    params.chamber_pressure  = chamber_pressure;
    params.grid_voltage      = grid_voltage;
    params.grid_current      = grid_current;

    // debug print params struct
    DEBUG("params.chamber_radius    = %0.3f meters  (%0.3f inches)\n", 
          params.chamber_radius,
          M_TO_IN(params.chamber_radius));
    DEBUG("params.grid_radius       = %0.3f meters  (%0.3f inches)\n", 
          params.grid_radius,
          M_TO_IN(params.grid_radius));
    DEBUG("params.chamber_pressure  = %0.3f pascal  (%0.1f mTorr)\n", 
          params.chamber_pressure,
          PASCAL_TO_MTORR(params.chamber_pressure));
    DEBUG("params.grid_voltage      = %0.0f volts\n", 
          params.grid_voltage);
    DEBUG("params.grid_current      = %0.3f amps\n", 
          params.grid_current);

    // loop over all cells
    // - init each cell
    // - update the volume of each shell
    // - keep track of the number of cells in the chamber
    int32_t x_idx, y_idx, z_idx, num_cell_in_chamber=0;
    float   q_grid = POTENTIAL_TO_CHARGE(params.grid_voltage, params.grid_radius);
    for (x_idx = 0; x_idx < MAX_CELL; x_idx++) {
        for (y_idx = 0; y_idx < MAX_CELL; y_idx++) {
            for (z_idx = 0; z_idx < MAX_CELL; z_idx++) {
                cell_t * c = &cell[x_idx][y_idx][z_idx];
                float    x, y, z, r;
                int32_t  shell_idx;

                // get the location of the center of the cell
                x = x_idx * CELL_SIZE - (MAX_CELL*CELL_SIZE/2) + CELL_SIZE/2;
                y = y_idx * CELL_SIZE - (MAX_CELL*CELL_SIZE/2) + CELL_SIZE/2;
                z = z_idx * CELL_SIZE - (MAX_CELL*CELL_SIZE/2) + CELL_SIZE/2;
                r = hypotenuse(x, y, z);

                // if this cell is not within the chamber then 
                // leave this cell uninitialized
                if (r >= params.chamber_radius) {
                    continue;
                }

                // init particle_list_head
                LIST_INIT(&c->particle_list_head);

                // init shell
                shell_idx = r / SHELL_SIZE;
                assert(shell_idx < MAX_SHELL);
                c->shell = &shell[shell_idx];

                // init the acceleration of the deuterium ion due to the 
                // electric field for this cell
                if (r > params.grid_radius) {
                    float f = COULOMB_FORCE(q_grid, PROTON_CHARGE, r);
                    c->ax = (f / D_MASS) * (x / r);
                    c->ay = (f / D_MASS) * (y / r);
                    c->az = (f / D_MASS) * (z / r);
                }

                // add the volume of the cell to the volume of the shell that contains it,
                // keep track of the num_cell_in_chamber
                c->shell->volume += CELL_VOLUME;
                num_cell_in_chamber++;

                // xxx
                if (shell_idx + 1 > max_shell) {
                    max_shell = shell_idx + 1;
                }
            }
        }
    }

    // init global variables for the kinetic velocity of D atoms at room temp
    roomtemp_d_velocity = TEMPERATURE_TO_VELOCITY(ROOM_TEMPERATURE_K, D_MASS);
    roomtemp_d_velocity_squared = roomtemp_d_velocity * roomtemp_d_velocity;

    // init particles
    int32_t i;
    int64_t start_us;
    max_particles = num_cell_in_chamber * AVG_SIM_PARTICLES_PER_CELL;
    assert(max_particles <= MAX_PARTICLES);
    DEBUG("initializing %d particles ...\n", max_particles);
    start_us = microsec_timer();
    for (i = 0; i < max_particles; i++) {
        init_particle(&particles[i]);
    }
    DEBUG("done initializing particles, %ld us\n", microsec_timer() - start_us);

    // determine num_real_particles_in_cell using the ideal gas law; and
    // determine num_real_particles_per_sim_particle
    num_real_particles_in_cell = 
        NUMBER_OF_PARTICLES(params.chamber_pressure, CELL_VOLUME, ROOM_TEMPERATURE_K);
    num_real_particles_per_sim_particle = num_real_particles_in_cell / AVG_SIM_PARTICLES_PER_CELL;

    // debug print 
    DEBUG("roomtemp_d_velocity = %f\n", roomtemp_d_velocity);
    DEBUG("num_cell_in_chamber = %d\n", num_cell_in_chamber);
    DEBUG("max_particles       = %d\n", max_particles);
    DEBUG("max_shell           = %d\n", max_shell);
    DEBUG("num_real_particles_in_cell          = %.0f  %e\n", 
          num_real_particles_in_cell, num_real_particles_in_cell);
    DEBUG("num_real_particles_per_sim_particle = %.0f  %e\n", 
          num_real_particles_per_sim_particle, num_real_particles_per_sim_particle);

    // init time
    time_model_secs = 0;

    // create the model_thread, and
    // wait for the model_thread to set its started flag
    pthread_t thread_id;
    pthread_create(&thread_id, NULL, model_thread, NULL);
    while (!model_thread_started) {
        usleep(1000);
    }

    // debug print init complete message
    INFO("---------- initialization complete ----------\n");
}

static void init_particle(particle_t * p)
{
    float x, y, z, vx, vy, vz;
    cell_t *c;

    // get a random location within the chamber
    random_location_within_chamber(&x, &y, &z, &c);

    // get a random velocity for this particle;
    // the direction is random, the magnitude is room temperature for D atom
    random_vector(roomtemp_d_velocity, &vx, &vy, &vz);

    // init particle position and velocity fields
    p->x = x;
    p->y = y;
    p->z = z;
    p->vx = vx;
    p->vy = vy;
    p->vz = vz;
    p->v = roomtemp_d_velocity;
    p->v_squared = roomtemp_d_velocity_squared;

    // init particle cell_entries field,
    // add the particle to the cell that it is contained in
    LIST_INSERT_HEAD(&c->particle_list_head, p, cell_entries);

    // init the remaining particle fields, 
    // except for the work_entries field which is left as null XXX maybe it could go on the queue here
    p->ion = false;
    p->cell = c;
    p->time_last_processed = 0;

    // updata shell info associated with this particle
    c->shell->number_of_atoms++;
}

// -----------------  MODEL CONTROL API  -------------------------------------------------

void model_start(void)
{
    // set run_request flag, and wait for ack
    run_request = true;
    while (!running) {
        usleep(1000);
    }
}

void model_stop(void)
{
    // clear run_request flag, and wait for ack
    run_request = false;
    while (running) {
        usleep(1000);
    }
}

bool model_is_running(void) 
{
    return running;
}

void model_terminate(void)
{
    model_stop();
}

// XXX much more work here
void model_get_data(pane_hndlr_display_graph_params_t ** graph_temp)
{
    pane_hndlr_display_graph_params_t *g;
    int32_t max_points_needed = 360;
    int32_t i;

    if (*graph_temp == NULL || (*graph_temp)->max_points_alloced < max_points_needed) {
        INFO("XXXXXXXXXXXXXXXX ALLOCING \n");
        *graph_temp = malloc(sizeof(pane_hndlr_display_graph_params_t) +
                             max_points_needed * sizeof(struct pane_hndlr_display_graph_point_s));
        g = *graph_temp;
        strcpy(g->title_str, "TITLE");
        strcpy(g->x_units_str, "X_UNIT");
        strcpy(g->y_units_str, "Y_UNIT");
        g->x_min = 0;
        g->x_max = 2 * M_PI;
        g->y_min = -1;
        g->y_max = 1;
        g->max_points_alloced = max_points_needed;
    }

    g = *graph_temp;
    for (i = 0; i < max_points_needed; i++) {
        g->points[i].x = i*(2*M_PI/360);
        g->points[i].y = sin(i*(2*M_PI/360));
    }
    g->max_points = max_points_needed;
}

// -----------------  MODEL_THREAD  ----------------------------------------------------------

#define MAX_WORK_LIST_HEAD 5000

static LIST_HEAD(wlh_s, particle_s) work_list_head[MAX_WORK_LIST_HEAD];
static int32_t wlh_idx;

static void schedule_work(particle_t *p, float delta_t);
static void process_particle(particle_t *p);

static void * model_thread(void * cx)
{
    struct wlh_s * wlh;
    particle_t * p;
    int32_t i;

    DEBUG("initializing\n");  // XXX move to init routine
    for (i = 0; i < MAX_WORK_LIST_HEAD; i++) {
        LIST_INIT(&work_list_head[i]);
    }
    for (i = 0; i < max_particles; i++) {
        LIST_INSERT_HEAD(&work_list_head[i%MAX_WORK_LIST_HEAD], &particles[i], work_entries);
    }

    DEBUG("starting\n");
    model_thread_started = true;

    while (true) {
        // wait for run_request to be set
        while (run_request == false) {
            running = false;
            usleep(1000);
        }
        running = true;

        // process the particles on this list;
        // when done all particles will have been moved to other work lists
        wlh = &work_list_head[wlh_idx];
        while ((p = wlh->lh_first) != NULL) {
            process_particle(p);
        }

        // increment the time, and the wlh_idx
        time_model_secs += MODEL_TIME_INCREMENT;
        wlh_idx = wlh_idx + 1;
        if (wlh_idx == MAX_WORK_LIST_HEAD) {
            wlh_idx = 0;
        }
    }

    DEBUG("terminating\n");
    return NULL;
}

static void schedule_work(particle_t *p, float delta_t)
{
    struct wlh_s * future_wlh;
    int32_t n;

    // remove p from the work list it is currently on
    LIST_REMOVE(p, work_entries);

    // add p to work list delta_t in the future
    n = delta_t / MODEL_TIME_INCREMENT;
    assert(n < MAX_WORK_LIST_HEAD);
    future_wlh = &work_list_head[(wlh_idx + n) % MAX_WORK_LIST_HEAD];
    LIST_INSERT_HEAD(future_wlh, p, work_entries);
}

static void process_particle(particle_t *p)
{
    float t, x, y, z, vx, vy, vz, nd, mfp, f, new_v, new_v_squared;
    cell_t *cur_cell, * new_cell;
    shell_t *cur_shell, *new_shell;
    particle_t *ptgt;

    // XXX try prefetch the particle, and the ptgt
    // xxx summary comment

    // determine the time since this particle's last interaction 
    t = time_model_secs - p->time_last_processed;

    // determine the new location of this particle
    x = p->x + p->vx * t;
    y = p->y + p->vy * t;
    z = p->z + p->vz * t;

    // xxx
    cur_cell  = p->cell;
    cur_shell = cur_cell->shell;
    new_cell  = get_cell(x, y, z);

    // if the new particle location (x,y,z) is not within the chamber then
    // - pick a new velocity vector that is within the chamber .025 m
    // - xxx  update the comments here
    // - return
    if (new_cell == NULL) {
        while (true) {
            random_vector(roomtemp_d_velocity, &vx, &vy, &vz);
            t = .020 / roomtemp_d_velocity;  // xxx tbd, this extends out of chamber
            x = p->x + vx * t;
            y = p->y + vy * t;
            z = p->z + vz * t;
            if (get_cell(x, y, z) != NULL) {
                p->vx = vx;
                p->vy = vy;
                p->vz = vz;
                p->v = roomtemp_d_velocity;
                p->v_squared = roomtemp_d_velocity_squared;
                p->time_last_processed = time_model_secs;
                nd = cur_shell->number_of_atoms * num_real_particles_per_sim_particle / cur_shell->volume;
                mfp = MEAN_FREE_PATH(H2_CROSS_SECTION,nd);
                schedule_work(p, mfp/p->v); 
                return;
            }
        }
    }

    // xxx
    new_shell = new_cell->shell;

    // choose a random particle in the new_cell;
    // ensure it is different than p
    // xxx if (__glibc_unlikely(new_cell->radius_idx >= max_shell)) {
    ptgt = new_cell->particle_list_head.lh_first;
    if (ptgt == p) {
        ptgt = LIST_NEXT(p, cell_entries);
    }

    // xxx improve this debug code
    static uint64_t okay, notokay;
    if (ptgt) {
        okay++;
    } else {
        notokay++;
    }
    if ((okay % 1000000) == 0) {
        DEBUG("XXX %12ld %12ld : %f\n", okay, notokay, (double)notokay/okay);
    }

    // both p and ptgt will be given the same velocity, by conservation of energy;
    // determine this new velocity by taking the average of the p and ptgt velocity squared
    if (ptgt) {
        new_v_squared = 0.5 * (p->v_squared + ptgt->v_squared);
        new_v = sqrtf(new_v_squared);
    } else {
        new_v = p->v;
        new_v_squared = p->v_squared;
    }

    // xxx optimize
    nd = new_shell->number_of_atoms * num_real_particles_per_sim_particle / new_shell->volume;
    mfp = MEAN_FREE_PATH(H2_CROSS_SECTION,nd);
    t = mfp / new_v;

    // update p
    p->x = x;
    p->y = y;
    p->z = z;
    f = new_v / p->v;
    p->vx = p->vx * f;
    p->vy = p->vy * f;
    p->vz = p->vz * f;
    p->v = new_v;
    p->v_squared = new_v_squared;
    p->time_last_processed = time_model_secs;
    if (new_cell != cur_cell) {
        LIST_REMOVE(p, cell_entries);
        LIST_INSERT_HEAD(&new_cell->particle_list_head, p, cell_entries);
        p->cell = new_cell;
    }
    schedule_work(p, t);

    // update ptgt
    if (ptgt) {
        f = new_v / ptgt->v;
        ptgt->vx = ptgt->vx * f;
        ptgt->vy = ptgt->vy * f;
        ptgt->vz = ptgt->vz * f;
        ptgt->v = new_v;
        ptgt->v_squared = new_v_squared;
        ptgt->time_last_processed = time_model_secs;
        schedule_work(ptgt, t);
    }
}

