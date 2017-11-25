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
#include <pthread.h>
#include <sys/mman.h>

#include "model.h"
#include "physics.h"
#include "util_misc.h"

//
// defines
//

// model time increment
#define DELTA_T_NS  1L
#define DELTA_T_SECS (1e-9 * DELTA_T_NS)

// 
// typedefs
//

//
// variables
//

static bool run_request;
static bool running;
 
//
// prototypes
//

void init_particle(particle_t * p);
void * model_thread(void * cx);

//
// inline functions
//

// xxx comment
inline cell_t * get_cell(int32_t x_nm, int32_t y_nm, int32_t z_nm)
{
    // xxx return null if cell is not in chamber ??
    int32_t x_idx = (x_nm + (MAX_CELL*CELL_SIZE_NM/2)) / CELL_SIZE_NM;
    int32_t y_idx = (y_nm + (MAX_CELL*CELL_SIZE_NM/2)) / CELL_SIZE_NM;
    int32_t z_idx = (z_nm + (MAX_CELL*CELL_SIZE_NM/2)) / CELL_SIZE_NM;
    assert(x_idx >= 0 && x_idx < MAX_CELL);
    assert(y_idx >= 0 && y_idx < MAX_CELL);
    assert(z_idx >= 0 && z_idx < MAX_CELL);
    
    return &cell[x_idx][y_idx][z_idx];
}

// xxx comment
inline void random_location_within_chamber(int32_t *x_nm_arg, int32_t *y_nm_arg, int32_t *z_nm_arg)
{
    int32_t x_nm, y_nm, z_nm;
    cell_t * c;

    while (true) {
        x_nm = random_range(-params.chamber_radius_nm, params.chamber_radius_nm);
        y_nm = random_range(-params.chamber_radius_nm, params.chamber_radius_nm);
        z_nm = random_range(-params.chamber_radius_nm, params.chamber_radius_nm);
        c = get_cell(x_nm, y_nm, z_nm);
        if (c->shell) {
            *x_nm_arg = x_nm;
            *y_nm_arg = y_nm;
            *z_nm_arg = z_nm;
            break;
        }
    }
}

// -----------------  MODEL INIT FROM PARAMS  --------------------------------------------

void model_init_from_params(char * params_str)
{
    int32_t ret, cnt, i;
    double chamber_diameter_mm, grid_diameter_mm, chamber_pressure_mtorr, grid_voltage_kv, grid_current_ma;
    int32_t x_idx, y_idx, z_idx;
    int32_t x_nm, y_nm, z_nm;
    double q_grid;
    int64_t num_cell_in_chamber;
    int64_t num_real_particles_in_cell;
    pthread_t thread_id;

    // debug print sizes
    DEBUG("MAX_PARTICLES     = %d\n", MAX_PARTICLES);
    DEBUG("MAX_CELL          = %d\n", MAX_CELL);
    DEBUG("MAX_SHELL         = %d\n", MAX_SHELL);
    DEBUG("sizeof(particles) = %ld MB\n", sizeof(particles) / MB);
    DEBUG("sizeof(cell)      = %ld MB\n", sizeof(cell) / MB);
    DEBUG("sizeof(shell)     = %ld MB\n", sizeof(shell) / MB);

    // don't dump the particles and cell arrays, because they are so large
    ret = madvise((void*)ROUND_UP(particles,PAGE_SIZE), sizeof(particles)-PAGE_SIZE, MADV_DONTDUMP);
    if (ret != 0) {
        FATAL("madvise particles ret %d, %s\n", ret, strerror(errno));
    }
    ret = madvise((void*)ROUND_UP(cell,PAGE_SIZE), sizeof(cell)-PAGE_SIZE, MADV_DONTDUMP);
    if (ret != 0) {
        FATAL("madvise cell ret %d, %s\n", ret, strerror(errno));
    }

    // convert param_str to params
    INFO("params_str = %s\n", params_str);
    cnt = sscanf(params_str, 
                 "%lf,%lf,%lf,%lf,%lf",
                 &chamber_diameter_mm,
                 &grid_diameter_mm,
                 &chamber_pressure_mtorr,
                 &grid_voltage_kv,
                 &grid_current_ma);
    if (cnt != 5) {
        FATAL("params_str '%s' must contain 5 values\n", params_str);
    }

    if (chamber_diameter_mm < MIN_CHAMBER_DIAMETER_MM ||
        chamber_diameter_mm > MAX_CHAMBER_DIAMETER_MM) {
        FATAL("chamber_diameter_mm %lf is out of range\n", chamber_diameter_mm);
    }
    if (grid_diameter_mm < MIN_GRID_DIAMETER_MM ||
        grid_diameter_mm > MAX_GRID_DIAMETER_MM(chamber_diameter_mm)) {
        FATAL("grid_diameter_mm %lf is out of range\n", grid_diameter_mm);
    }
    if (chamber_pressure_mtorr < MIN_CHAMBER_PRESSURE_MTORR ||
        chamber_pressure_mtorr > MAX_CHAMBER_PRESSURE_MTORR) {
        FATAL("chamber_pressure_mtorr %lf is out of range\n", chamber_pressure_mtorr);
    }
    if (grid_voltage_kv > 0 || 
        grid_voltage_kv < MAX_GRID_VOLTAGE_KV)  {
        FATAL("grid_voltage_kv %lf is out of range\n", grid_voltage_kv);
    }
    if (grid_current_ma < 0 ||
        grid_current_ma > MAX_GRID_CURRENT_MA) {
        FATAL("grid_current_ma %lf is out of range\n", grid_current_ma);
    }

    params.chamber_radius_nm      = MM_TO_NM(chamber_diameter_mm / 2 / SHELL_SIZE_MM * SHELL_SIZE_MM);
    params.grid_radius_nm         = MM_TO_NM(grid_diameter_mm / 2);
    params.chamber_pressure_utorr = chamber_pressure_mtorr * 1000;
    params.grid_voltage_v         = grid_voltage_kv * 1000;
    params.grid_current_ua        = grid_current_ma * 1000;

    // debug print params
    DEBUG("params.chamber_radius_nm      = %d (%f inches)\n", 
          params.chamber_radius_nm,
          NM_TO_IN(params.chamber_radius_nm));
    DEBUG("params.grid_radius_nm         = %d (%f inches)\n", 
          params.grid_radius_nm,
          NM_TO_IN(params.grid_radius_nm));
    DEBUG("params.chamber_pressure_utorr = %d\n", params.chamber_pressure_utorr);
    DEBUG("params.grid_voltage_v         = %d\n", params.grid_voltage_v);
    DEBUG("params.grid_current_ua        = %d\n", params.grid_current_ua);

    // init max_shell 
    max_shell = params.chamber_radius_nm / SHELL_SIZE_NM;
    DEBUG("max_shell = %d\n", max_shell);
    assert(max_shell <= MAX_SHELL);

    // loop over all cells
    // - init each cell
    // - update the volume of each shell
    // - keep track of the number of cells in the chamber
    num_cell_in_chamber = 0;
    q_grid = POTENTIAL_TO_CHARGE(params.grid_voltage_v, NM_TO_M((double)params.grid_radius_nm));
    for (x_idx = 0; x_idx < MAX_CELL; x_idx++) {
        for (y_idx = 0; y_idx < MAX_CELL; y_idx++) {
            for (z_idx = 0; z_idx < MAX_CELL; z_idx++) {
                cell_t * c = &cell[x_idx][y_idx][z_idx];

                // get the location of the center of the cell
                x_nm = x_idx * CELL_SIZE_NM - (MAX_CELL*CELL_SIZE_NM/2) + CELL_SIZE_NM/2;
                y_nm = y_idx * CELL_SIZE_NM - (MAX_CELL*CELL_SIZE_NM/2) + CELL_SIZE_NM/2;
                z_nm = z_idx * CELL_SIZE_NM - (MAX_CELL*CELL_SIZE_NM/2) + CELL_SIZE_NM/2;

                // init particle_list_head
                LIST_INIT(&c->particle_list_head);

                // init r_nm (distance from center of chamber to center of this cell)
                c->r_nm = hypotenuse(x_nm, y_nm, z_nm);

                // init shell, this is set to NULL if the cell is outside the chamber
                if (c->r_nm < params.chamber_radius_nm) {
                    int32_t idx = c->r_nm / SHELL_SIZE_NM;
                    assert(idx < max_shell);
                    c->shell = &shell[idx];
                } else {
                    c->shell = NULL;
                }

                // init xa_nmperns2 (and y,z); this is the acceleration of the
                // deuterium ion due to the electric field for this cell
                if (c->r_nm > params.grid_radius_nm) {
                    double f = COULOMB_FORCE(q_grid, PROTON_CHARGE, NM_TO_M((double)c->r_nm));
                    double m = AMU_TO_KG(D_AMU);
                    c->xa_nmperns2 = (f / m) * ((double)x_nm / c->r_nm) * 1e-9;
                    c->ya_nmperns2 = (f / m) * ((double)y_nm / c->r_nm) * 1e-9;
                    c->za_nmperns2 = (f / m) * ((double)z_nm / c->r_nm) * 1e-9;
                } else {
                    c->xa_nmperns2 = 0;
                    c->ya_nmperns2 = 0;
                    c->za_nmperns2 = 0;
                }

                // if this cell is within the chamber then
                //   add the volume of the cell to the volume of the shell that contains it,
                //   keep track of the num_cell_in_chamber
                // endif
                if (c->shell != NULL) {
                    c->shell->volume_cu_mm += CELL_VOLUME_CU_MM;
                    num_cell_in_chamber++;
                }
            }
        }
    }
    DEBUG("num_cell_in_chamber = %ld\n", num_cell_in_chamber);

    // init max_particles
    max_particles = num_cell_in_chamber * AVG_SIM_PARTICLES_PER_CELL;
    DEBUG("max_particles = %d\n", max_particles);
    assert(max_particles <= MAX_PARTICLES);

    // init particles
    DEBUG("initializing particles ...\n");
    int64_t start_us = microsec_timer();
    for (i = 0; i < max_particles; i++) {
        init_particle(&particles[i]);
    }
    DEBUG("done initializing particles, %ld us\n", microsec_timer() - start_us);

    // determine num_real_particles_per_sim_particle 
    // using the ideal gas law
    num_real_particles_in_cell = 
        NUMBER_OF_PARTICLES(UTORR_TO_PASCAL(params.chamber_pressure_utorr),
                            CELL_VOLUME_CU_MM/1e9,
                            ROOM_TEMPERATURE_K);
    num_real_particles_per_sim_particle = (double)num_real_particles_in_cell / AVG_SIM_PARTICLES_PER_CELL;
    DEBUG("num_real_particles_in_cell          = %ld  %e\n", 
          num_real_particles_in_cell, (double)num_real_particles_in_cell);
    DEBUG("num_real_particles_per_sim_particle = %lf  %e\n", 
          num_real_particles_per_sim_particle, num_real_particles_per_sim_particle);

    // init time_ns
    time_ns = 0;

    // create the model_thread
    pthread_create(&thread_id, NULL, model_thread, NULL);
}

void init_particle(particle_t * p)
{
    int32_t x_nm, y_nm, z_nm;
    int32_t xv_nmperns, yv_nmperns, zv_nmperns;
    cell_t *c;

    static int32_t roomtemp_d_velocity_mpers;
    static int32_t roomtemp_d_velocity_nmperns;

    // on first call initialize roomtemp_d_velocity_mpers
    if (roomtemp_d_velocity_mpers == 0) {
        roomtemp_d_velocity_mpers = TEMPERATURE_TO_VELOCITY(ROOM_TEMPERATURE_K, D_MASS);
        roomtemp_d_velocity_nmperns = MPERS_TO_NMPERNS(roomtemp_d_velocity_mpers);
        DEBUG("roomtemp_d_velocity_mpers = %d\n", roomtemp_d_velocity_mpers);
    }

    // get a random location within the chamber
    random_location_within_chamber(&x_nm, &y_nm, &z_nm);

    // get a random velocity for this particle;
    // the direction is random, the magnitude is room temperature for D atom
    random_vector(roomtemp_d_velocity_nmperns, &xv_nmperns, &yv_nmperns, &zv_nmperns);

    // init particle fields
    p->x_nm = x_nm;
    p->y_nm = y_nm;
    p->z_nm = z_nm;
    p->xv_nmperns = xv_nmperns;
    p->yv_nmperns = yv_nmperns;
    p->zv_nmperns = zv_nmperns;
    p->velocity_squared_m2pers2 = (int64_t)roomtemp_d_velocity_mpers * (int64_t)roomtemp_d_velocity_mpers;
    p->ion = false;

    // add the particle to the cell that it is contained in
    c = get_cell(x_nm, y_nm, z_nm);
    LIST_INSERT_HEAD(&c->particle_list_head, p, entries);

    // updata shell info associated with this particle
    c->shell->number_of_atoms++;
    c->shell->sum_velocity_squared_m2pers2 +=  p->velocity_squared_m2pers2;
}

// -----------------  MODEL INIT FROM FILE  ----------------------------------------------

void model_init_from_file(char * filename_str)
{
    // xxx tbd
}

// -----------------  RUN CONTROL API  ---------------------------------------------------

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

void model_terminate(void)
{
    model_stop();
}

// -----------------  MODEL THREADS  ---------------------------------------------------------

void * model_thread(void * cx)
{
    DEBUG("starting\n");
    while (true) {
        // wait for run_request to be set
        while (run_request == false) {
            running = false;
            usleep(1000);
        }
        running = true;

        sleep(1);
#if 0
        // get the xxx list head

        // process the particles on this list;
        // when done all particles will have been moved to other list
        for ...
            process_particle


        // this list has been processed, so clear the list head


        // increment the time
        time_ns += DELTA_T_NS;
#endif

    }
    DEBUG("terminating\n");
    return NULL;
}

#if 0
void process_particle(particle_t * p)
{
    if (p->ion) {
        p->xv_nmperns += cell->dxv_nmperns;
        p->yv_nmperns += cell->dyv_nmperns;
        p->zv_nmperns += cell->dzv_nmperns;
    }

    // update the particle's position, and 
    // determine the cell it is now in following the position update
    p->x_nm += p->xv_nmperns;
    p->y_nm += p->yv_nmperns;
    p->z_nm += p->zv_nmperns;
    new_cell = get_cell(p->x_nm, p->y_nm, p->z_nm);

    // if the particle is now outside the chamber then 
    // - put the particle back where it was, and
    // - update its velocity to be roomtemp toward the center of the chamber
    if (__glibc_unlikely(new_cell->radius_idx >= max_shell)) {
    }

    // if the particle has moved to another cell then
    // remove it from this cell list, and add it to the
    // new cell list
    if (__glibc_unlikely(new_cell != cell)) {
        LIST_REMOVE(p, entries);
        LIST_INSERT_HEAD(&cell->particle_list_head, p, entries);

        // if the new cell is contained in a different radius shell then
        // update the old and new radius shells
        if (new_cell->radius_idx != cell->radius_idx) {
            radius_t * r_old = &radius[cell->radius_idx];
            if (!p->ion) {
                r_old->number_of_atoms--;
            } else {
                r_old->number_of_ions--;
            }
            r_old->sum_velocity_squared_m2pers2 -= p->velocity_squared_m2pers2;

            radius_t * r_new = &radius[new_cell->radius_idx];
            if (!p->ion) {
                r_new->number_of_atoms++;
            } else {
                r_new->number_of_ions++;
            }
            r_new->sum_velocity_squared_m2pers2 += p->velocity_squared_m2pers2;

            stats.moved_radius++;  // xxx is this needed
        }

    } 
}
#endif

#if 0
static struct {
    int64_t start_us;
    int64_t done_us;
    int64_t processed;
    int64_t moved_cell;
    int64_t moved_radius;
} stats;

void * control_thread(void * cx)
{
    while (true) {
        // wait for run_request to be set
        if (!run_request) {
            running = false;
            while (run_request == false) {
                usleep(1000);
            }
            running = true;
        }

        // initialize for this cycle
        time_ns += DELTA_T_NS;
        cell_list_idx = 0;
        memset(&stats,0,sizeof(stats));

        // release the work threads
        DEBUG("cycle start: time_ns = %ld\n", time_ns);
        stats.start_us = microsec_timer();
        pthread_barrier_wait(&barrier);

        // wait for the work threads to complete
        pthread_barrier_wait(&barrier);
        stats.done_us = microsec_timer();
        DEBUG("cycle complete: duration = %ld  processed = %ld  moved_cell = %ld  moved_radius = %ld\n",
              stats.done_us - stats.start_us,
              stats.processed,
              stats.moved_cell,
              stats.moved_radius);
    }
}

void * work_thread(void * cx)
{
    int32_t      thread_id  __attribute__ ((unused)) = (int64_t)cx;
    int32_t      idx;
    cell_t   * cell;
    cell_t   * new_cell;
    particle_t * p;
    particle_t * p_next;
    LIST_HEAD(head_s, particle_s) moved_list_head;

    DEBUG("starting thread_id = %d\n", thread_id);

    while (true) {
        // wait to be started 
        pthread_barrier_wait(&barrier);

        // loop until all cell have been processed
        // - get the next location box 
        // - process all the particles within the location box
        //   that have not already been processed 
        // xxx try prefetch to improve performance
        while (true) {
            // get the next location box, if no more then break
            idx = __sync_fetch_and_add(&cell_list_idx, 1);
            if (idx >= max_cell_list) {
                break;
            }
            cell = cell_list[idx];

            // init moved_list_head; this is used to save a list of
            // the particles that will be moved from this cell, so that
            // these can be processed after releasing the particle_list_spinlock
            LIST_INIT(&moved_list_head);

            // acquire spinlock to protect access to the cell's particle_list
            pthread_spin_lock(&cell->particle_list_spinlock);
    
            // loop over all particles within the location box,
            // and process them if they have not already been
            // processed; a particle may already have been processed
            // if it has just been moved into this location box
            for (p = cell->particle_list_head.lh_first; p != NULL; p = p_next) {
                // save ptr to the next particle in the list
                p_next = p->entries.le_next;

                // if the particle has already been processed then continue
                if (p->last_processed_time_ns == time_ns) {
                    continue;
                }

x               // xxx velocity
x               if (p->ion) {
x                   p->xv_nmperns += cell->dxv_nmperns;
x                   p->yv_nmperns += cell->dyv_nmperns;
x                   p->zv_nmperns += cell->dzv_nmperns;

-                   int32_t v_nmperns;  // xxx DEBUG ....
-                   static int32_t v_max_nmperns;
-                   v_nmperns = hypotenuse( p->xv_nmperns, p->yv_nmperns, p->zv_nmperns);
-                   if (v_nmperns > v_max_nmperns) {
-                       v_max_nmperns = v_nmperns;
-                       //DEBUG("xxx v_max_nmperns = %d\n", v_max_nmperns);
-                   }
-               }

x               // update the particle's position, and 
x               // determine the cell it is now in following the position update
x               p->x_nm += p->xv_nmperns;
x               p->y_nm += p->yv_nmperns;
x               p->z_nm += p->zv_nmperns;
x               new_cell = get_cell(p->x_nm, p->y_nm, p->z_nm);

x               // if the particle is now outside the chamber then 
x               // - put the particle back where it was, and
x               // - update its velocity to be roomtemp toward the center of the chamber
x               if (__glibc_unlikely(new_cell->radius_idx >= max_shell)) {
                    // put the particle back where it was
                    p->x_nm -= p->xv_nmperns;
                    p->y_nm -= p->yv_nmperns;
                    p->z_nm -= p->zv_nmperns;
                    new_cell = cell;

                    // update its velocity to be roomtemp toward the center of the chamber
                    p->xv_nmperns = (int64_t)(-p->x_nm) * roomtemp_d_velocity_nmperns / cell->r_nm;
                    p->yv_nmperns = (int64_t)(-p->y_nm) * roomtemp_d_velocity_nmperns / cell->r_nm;
                    p->zv_nmperns = (int64_t)(-p->z_nm) * roomtemp_d_velocity_nmperns / cell->r_nm;

#if 0 // xxx debug, delete
                    static bool first = true;
                    if (first) {
                        first = false;
                        int32_t vel = hypotenuse(p->xv_nmperns, p->yv_nmperns, p->zv_nmperns);
                        DEBUG("vel = %d  roomtemp_d_vel = %d   vel_vect= %d %d %d  pos = %d %d %d\n",
                              vel, roomtemp_d_velocity_nmperns, 
                              p->xv_nmperns, p->yv_nmperns, p->zv_nmperns,
                              p->x_nm, p->y_nm, p->z_nm);
                    }
#endif
                }

x               // if the particle has moved to another cell then
x               // remove it from this cell list, and add it to the
x               // moved_list; the reason for the moved_list is that we
x               // can not put this particle on the new_cell particle_list 
x               // at this time due to a possible deadlock; we first 
x               // need to drop the particle_list_spinlock
x               if (__glibc_unlikely(new_cell != cell)) {
x                   LIST_REMOVE(p, entries);
x                   LIST_INSERT_HEAD(&moved_list_head, p, entries);
x               } 

                // indicate that this particle has been processed
                p->last_processed_time_ns = time_ns;
                __sync_fetch_and_add(&stats.processed,1);
                    
            }

            // unlock the particle_list_spinlock
            pthread_spin_unlock(&cell->particle_list_spinlock);

            // process the moved_list
            while ((p = moved_list_head.lh_first) != NULL) {
                LIST_REMOVE(p, entries);

                // determine the new cell that the particle is being moved to
                new_cell = get_cell(p->x_nm, p->y_nm, p->z_nm);

x              // if the new cell is contained in a different radius shell then
x               // update the old and new radius shells
x               if (new_cell->radius_idx != cell->radius_idx) {
x                   radius_t * r_old = &radius[cell->radius_idx];
x                   if (!p->ion) {
x                       __sync_fetch_and_sub(&r_old->number_of_atoms,1);
x                   } else {
x                       __sync_fetch_and_sub(&r_old->number_of_ions,1);
x                   }
x                   __sync_fetch_and_sub(&r_old->sum_velocity_squared_m2pers2, p->velocity_squared_m2pers2);

x                   radius_t * r_new = &radius[new_cell->radius_idx];
x                   if (!p->ion) {
x                       __sync_fetch_and_add(&r_new->number_of_atoms,1);
x                   } else {
x                       __sync_fetch_and_add(&r_new->number_of_ions,1);
x                   }
x                   __sync_fetch_and_add(&r_new->sum_velocity_squared_m2pers2, p->velocity_squared_m2pers2);

x                   __sync_fetch_and_add(&stats.moved_radius,1);
x               }

                // add the particle to the particle_list in the new cell
                pthread_spin_lock(&new_cell->particle_list_spinlock);
                LIST_INSERT_HEAD(&new_cell->particle_list_head, p, entries);
                pthread_spin_unlock(&new_cell->particle_list_spinlock);

                // stats update
                __sync_fetch_and_add(&stats.moved_cell,1);
            }
        }

        // indicate completed to the control_thread
        pthread_barrier_wait(&barrier);
    }
}
#endif
