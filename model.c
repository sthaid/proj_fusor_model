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

static bool    run_request;
static bool    running;
 
//
// prototypes
//

void init_particle(particle_t * p);
void * model_thread(void * cx);

//
// inline functions
//

inline locbox_t * get_locbox(int32_t x_nm, int32_t y_nm, int32_t z_nm)
{
    int32_t x_idx = (x_nm + (MAX_LOCBOX*LOCBOX_SIZE_NM/2)) / LOCBOX_SIZE_NM;
    int32_t y_idx = (y_nm + (MAX_LOCBOX*LOCBOX_SIZE_NM/2)) / LOCBOX_SIZE_NM;
    int32_t z_idx = (z_nm + (MAX_LOCBOX*LOCBOX_SIZE_NM/2)) / LOCBOX_SIZE_NM;
    assert(x_idx >= 0 && x_idx < MAX_LOCBOX);
    assert(y_idx >= 0 && y_idx < MAX_LOCBOX);
    assert(z_idx >= 0 && z_idx < MAX_LOCBOX);
    
    return &locbox[x_idx][y_idx][z_idx];
}

inline int32_t hypotenuse(int32_t x, int32_t y, int32_t z)
{
    int64_t hypotenuse_squared = (int64_t)x * (int64_t)x +
                                 (int64_t)y * (int64_t)y +
                                 (int64_t)z * (int64_t)z;
    return sqrt(hypotenuse_squared);
}

// -----------------  MODEL INIT FROM PARAMS  --------------------------------------------

void model_init_from_params(char * params_str)
{
    int32_t ret, cnt, i;
    double chamber_diameter_mm, grid_diameter_mm, chamber_pressure_mtorr, grid_voltage_kv, grid_current_ma;
    int32_t x_idx, y_idx, z_idx;
    int32_t x_nm, y_nm, z_nm;
    double q_grid;
    int64_t num_locbox_in_chamber;
    int64_t num_real_particles_in_locbox;
    pthread_t thread_id;
    int64_t start_us;

    // don't dump the particles and locbox arrays, because they are so large
    ret = madvise((void*)ROUND_UP(particles,PAGE_SIZE), sizeof(particles)-PAGE_SIZE, MADV_DONTDUMP);
    if (ret != 0) {
        FATAL("madvise particles ret %d, %s\n", ret, strerror(errno));
    }
    ret = madvise((void*)ROUND_UP(locbox,PAGE_SIZE), sizeof(locbox)-PAGE_SIZE, MADV_DONTDUMP);
    if (ret != 0) {
        FATAL("madvise locbox ret %d, %s\n", ret, strerror(errno));
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

    params.chamber_radius_nm      = 
            (int32_t)chamber_diameter_mm / 2 / RADIUS_SHELL_SIZE_MM * RADIUS_SHELL_SIZE_MM * 1000000;
    params.grid_radius_nm         = grid_diameter_mm * 1000000 / 2;
    params.chamber_pressure_utorr = chamber_pressure_mtorr * 1000;
    params.grid_voltage_v         = grid_voltage_kv * 1000;
    params.grid_current_ua        = grid_current_ma * 1000;

    DEBUG("params.chamber_radius_nm      = %d (%f inches)\n", 
          params.chamber_radius_nm,
          M_TO_IN(params.chamber_radius_nm / 1e9));
    DEBUG("params.grid_radius_nm         = %d (%f inches)\n", 
          params.grid_radius_nm,
          M_TO_IN(params.grid_radius_nm / 1e9));
    DEBUG("params.chamber_pressure_utorr = %d\n", params.chamber_pressure_utorr);
    DEBUG("params.grid_voltage_v         = %d\n", params.grid_voltage_v);
    DEBUG("params.grid_current_ua        = %d\n", params.grid_current_ua);

    max_radius = params.chamber_radius_nm / (RADIUS_SHELL_SIZE_MM * 1000000L);
    DEBUG("max_radius = %d\n", max_radius);

    // init locbox , xxx call cell
    num_locbox_in_chamber = 0;
    q_grid = POTENTIAL_TO_CHARGE(params.grid_voltage_v, NM_TO_M((double)params.grid_radius_nm));
    for (x_idx = 0; x_idx < MAX_LOCBOX; x_idx++) {
        for (y_idx = 0; y_idx < MAX_LOCBOX; y_idx++) {
            for (z_idx = 0; z_idx < MAX_LOCBOX; z_idx++) {
                locbox_t * lb = &locbox[x_idx][y_idx][z_idx];

                // get the location of the center of the locbox
                x_nm = x_idx * LOCBOX_SIZE_NM - (MAX_LOCBOX*LOCBOX_SIZE_NM/2) + LOCBOX_SIZE_NM/2;
                y_nm = y_idx * LOCBOX_SIZE_NM - (MAX_LOCBOX*LOCBOX_SIZE_NM/2) + LOCBOX_SIZE_NM/2;
                z_nm = z_idx * LOCBOX_SIZE_NM - (MAX_LOCBOX*LOCBOX_SIZE_NM/2) + LOCBOX_SIZE_NM/2;

                // init particle_list_head
                LIST_INIT(&lb->particle_list_head);

                // init r_nm (the radius of the locbox center), and radius_idx
                lb->r_nm = hypotenuse(x_nm, y_nm, z_nm);
                lb->radius_idx = lb->r_nm / (RADIUS_SHELL_SIZE_MM * 1000000L);  // xxx macro
                if (lb->radius_idx >= MAX_RADIUS) {
                    FATAL("radius_idx %d too big MAX_RADIUS = %ld\n", lb->radius_idx, MAX_RADIUS);
                }

                if (lb->r_nm > params.grid_radius_nm) {
                    double f = COULOMB_FORCE(q_grid, PROTON_CHARGE, NM_TO_M((double)lb->r_nm));
                    double m = AMU_TO_KG(D_AMU);
                    lb->xa_nmperns2 = (f / m) * ((double)x_nm / lb->r_nm) * 1e-9;
                    lb->ya_nmperns2 = (f / m) * ((double)y_nm / lb->r_nm) * 1e-9;
                    lb->za_nmperns2 = (f / m) * ((double)z_nm / lb->r_nm) * 1e-9;
                } else {
                    lb->xa_nmperns2 = 0;
                    lb->ya_nmperns2 = 0;
                    lb->za_nmperns2 = 0;
                }

                // xxx comment
                radius[lb->radius_idx].volume_cu_mm += LOCBOX_VOLUME_CU_MM;
                if (lb->radius_idx < max_radius) {
                    num_locbox_in_chamber++;
                }
            }
        }
    }
    DEBUG("num_locbox_in_chamber = %ld\n", num_locbox_in_chamber);

    // init particles and max_particles
    max_particles = num_locbox_in_chamber * AVERAGE_PARTICLES_PER_LOCBOX;
    DEBUG("max_particles = %d\n", max_particles);
    DEBUG("initializing particles ...\n");
    start_us = microsec_timer();
    for (i = 0; i < max_particles; i++) {
        init_particle(&particles[i]);
    }
    DEBUG("done initializing particles, %ld us\n", microsec_timer() - start_us);

    // determine num_real_particles_per_virtual_particle 
    // using the ideal gas law
    num_real_particles_in_locbox = 
        NUMBER_OF_PARTICLES(UTORR_TO_PASCAL(params.chamber_pressure_utorr),
                            LOCBOX_VOLUME_CU_MM/1000000000.0,
                            ROOM_TEMPERATURE_K);
    num_real_particles_per_virtual_particle = num_real_particles_in_locbox / AVERAGE_PARTICLES_PER_LOCBOX;
    DEBUG("num_real_particles_in_locbox            = %ld  %e\n", 
          num_real_particles_in_locbox, (double)num_real_particles_in_locbox);
    DEBUG("num_real_particles_per_virtual_particle = %ld  %e\n", 
          num_real_particles_per_virtual_particle, (double)num_real_particles_per_virtual_particle);

    // init time_ns
    time_ns = 0;

    // create threads
    pthread_create(&thread_id, NULL, model_thread, NULL);
}

void init_particle(particle_t * p)
{
    int32_t x_nm, y_nm, z_nm;
    int32_t x_dir, y_dir, z_dir, xv_nmperns, yv_nmperns, zv_nmperns;
    locbox_t *lb;

    static int32_t roomtemp_d_velocity_mpers;
    static int32_t roomtemp_d_velocity_nmperns;

    if (roomtemp_d_velocity_mpers == 0) {
        roomtemp_d_velocity_mpers = TEMPERATURE_TO_VELOCITY(ROOM_TEMPERATURE_K, D_MASS);
        roomtemp_d_velocity_nmperns = MPERS_TO_NMPERNS(roomtemp_d_velocity_mpers);
        DEBUG("roomtemp_d_velocity_mpers = %d\n", roomtemp_d_velocity_mpers);
    }

    // get a random location within the chamber
    while (true) {
        x_nm = random_range(-params.chamber_radius_nm, params.chamber_radius_nm);
        y_nm = random_range(-params.chamber_radius_nm, params.chamber_radius_nm);
        z_nm = random_range(-params.chamber_radius_nm, params.chamber_radius_nm);
        lb = get_locbox(x_nm, y_nm, z_nm);
        if (lb->radius_idx < max_radius) {
            break;
        }
    }

    // get a random velocity 
    // - the direction is random, the magnitude is room temperature for D atom
    // - the random direction is determined by picking a random point within a sphere
    int32_t hypot;
    while (true) {
        x_dir = random_range(-100000,100000);
        y_dir = random_range(-100000,100000);
        z_dir = random_range(-100000,100000);
        hypot = hypotenuse(x_dir,y_dir,z_dir);
        if (hypot > 1000 && hypot < 100000) {
            break;
        }
    }
    xv_nmperns = x_dir * roomtemp_d_velocity_nmperns / hypot;
    yv_nmperns = y_dir * roomtemp_d_velocity_nmperns / hypot;
    zv_nmperns = z_dir * roomtemp_d_velocity_nmperns / hypot;

#if 1  
    // sanity check the random velocity
    int32_t check_velocity_nmperns = hypotenuse(xv_nmperns, yv_nmperns, zv_nmperns);
    if (abs(check_velocity_nmperns-roomtemp_d_velocity_nmperns) > 2) {
        FATAL("roomtemp_d_velocity_nmperns = %d  check_velocity_nmperns = %d\n", 
              roomtemp_d_velocity_nmperns, check_velocity_nmperns);
    }
#endif

    // init particle fields
    p->x_nm = x_nm;
    p->y_nm = y_nm;
    p->z_nm = z_nm;
    p->xv_nmperns = xv_nmperns;
    p->yv_nmperns = yv_nmperns;
    p->zv_nmperns = zv_nmperns;
    p->velocity_squared_m2pers2 = (int64_t)roomtemp_d_velocity_mpers * (int64_t)roomtemp_d_velocity_mpers;
    p->ion = false;

    // add the particle to the locbox that it is contained in
    LIST_INSERT_HEAD(&lb->particle_list_head, p, entries);

    // updata associated radius info
    assert(lb->radius_idx < max_radius);
    radius_t * r = &radius[lb->radius_idx];
    r->number_of_atoms++;
    r->sum_velocity_squared_m2pers2 +=  p->velocity_squared_m2pers2;
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
        p->xv_nmperns += lb->dxv_nmperns;
        p->yv_nmperns += lb->dyv_nmperns;
        p->zv_nmperns += lb->dzv_nmperns;
    }

    // update the particle's position, and 
    // determine the locbox it is now in following the position update
    p->x_nm += p->xv_nmperns;
    p->y_nm += p->yv_nmperns;
    p->z_nm += p->zv_nmperns;
    new_lb = get_locbox(p->x_nm, p->y_nm, p->z_nm);

    // if the particle is now outside the chamber then 
    // - put the particle back where it was, and
    // - update its velocity to be roomtemp toward the center of the chamber
    if (__glibc_unlikely(new_lb->radius_idx >= max_radius)) {
    }

    // if the particle has moved to another locbox then
    // remove it from this locbox list, and add it to the
    // new locbox list
    if (__glibc_unlikely(new_lb != lb)) {
        LIST_REMOVE(p, entries);
        LIST_INSERT_HEAD(&lb->particle_list_head, p, entries);

        // if the new locbox is contained in a different radius shell then
        // update the old and new radius shells
        if (new_lb->radius_idx != lb->radius_idx) {
            radius_t * r_old = &radius[lb->radius_idx];
            if (!p->ion) {
                r_old->number_of_atoms--;
            } else {
                r_old->number_of_ions--;
            }
            r_old->sum_velocity_squared_m2pers2 -= p->velocity_squared_m2pers2;

            radius_t * r_new = &radius[new_lb->radius_idx];
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
    int64_t moved_locbox;
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
        locbox_list_idx = 0;
        memset(&stats,0,sizeof(stats));

        // release the work threads
        DEBUG("cycle start: time_ns = %ld\n", time_ns);
        stats.start_us = microsec_timer();
        pthread_barrier_wait(&barrier);

        // wait for the work threads to complete
        pthread_barrier_wait(&barrier);
        stats.done_us = microsec_timer();
        DEBUG("cycle complete: duration = %ld  processed = %ld  moved_locbox = %ld  moved_radius = %ld\n",
              stats.done_us - stats.start_us,
              stats.processed,
              stats.moved_locbox,
              stats.moved_radius);
    }
}

void * work_thread(void * cx)
{
    int32_t      thread_id  __attribute__ ((unused)) = (int64_t)cx;
    int32_t      idx;
    locbox_t   * lb;
    locbox_t   * new_lb;
    particle_t * p;
    particle_t * p_next;
    LIST_HEAD(head_s, particle_s) moved_list_head;

    DEBUG("starting thread_id = %d\n", thread_id);

    while (true) {
        // wait to be started 
        pthread_barrier_wait(&barrier);

        // loop until all locbox have been processed
        // - get the next location box 
        // - process all the particles within the location box
        //   that have not already been processed 
        // xxx try prefetch to improve performance
        while (true) {
            // get the next location box, if no more then break
            idx = __sync_fetch_and_add(&locbox_list_idx, 1);
            if (idx >= max_locbox_list) {
                break;
            }
            lb = locbox_list[idx];

            // init moved_list_head; this is used to save a list of
            // the particles that will be moved from this locbox, so that
            // these can be processed after releasing the particle_list_spinlock
            LIST_INIT(&moved_list_head);

            // acquire spinlock to protect access to the locbox's particle_list
            pthread_spin_lock(&lb->particle_list_spinlock);
    
            // loop over all particles within the location box,
            // and process them if they have not already been
            // processed; a particle may already have been processed
            // if it has just been moved into this location box
            for (p = lb->particle_list_head.lh_first; p != NULL; p = p_next) {
                // save ptr to the next particle in the list
                p_next = p->entries.le_next;

                // if the particle has already been processed then continue
                if (p->last_processed_time_ns == time_ns) {
                    continue;
                }

x               // xxx velocity
x               if (p->ion) {
x                   p->xv_nmperns += lb->dxv_nmperns;
x                   p->yv_nmperns += lb->dyv_nmperns;
x                   p->zv_nmperns += lb->dzv_nmperns;

-                   int32_t v_nmperns;  // xxx DEBUG ....
-                   static int32_t v_max_nmperns;
-                   v_nmperns = hypotenuse( p->xv_nmperns, p->yv_nmperns, p->zv_nmperns);
-                   if (v_nmperns > v_max_nmperns) {
-                       v_max_nmperns = v_nmperns;
-                       //DEBUG("xxx v_max_nmperns = %d\n", v_max_nmperns);
-                   }
-               }

x               // update the particle's position, and 
x               // determine the locbox it is now in following the position update
x               p->x_nm += p->xv_nmperns;
x               p->y_nm += p->yv_nmperns;
x               p->z_nm += p->zv_nmperns;
x               new_lb = get_locbox(p->x_nm, p->y_nm, p->z_nm);

x               // if the particle is now outside the chamber then 
x               // - put the particle back where it was, and
x               // - update its velocity to be roomtemp toward the center of the chamber
x               if (__glibc_unlikely(new_lb->radius_idx >= max_radius)) {
                    // put the particle back where it was
                    p->x_nm -= p->xv_nmperns;
                    p->y_nm -= p->yv_nmperns;
                    p->z_nm -= p->zv_nmperns;
                    new_lb = lb;

                    // update its velocity to be roomtemp toward the center of the chamber
                    p->xv_nmperns = (int64_t)(-p->x_nm) * roomtemp_d_velocity_nmperns / lb->r_nm;
                    p->yv_nmperns = (int64_t)(-p->y_nm) * roomtemp_d_velocity_nmperns / lb->r_nm;
                    p->zv_nmperns = (int64_t)(-p->z_nm) * roomtemp_d_velocity_nmperns / lb->r_nm;

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

x               // if the particle has moved to another locbox then
x               // remove it from this locbox list, and add it to the
x               // moved_list; the reason for the moved_list is that we
x               // can not put this particle on the new_lb particle_list 
x               // at this time due to a possible deadlock; we first 
x               // need to drop the particle_list_spinlock
x               if (__glibc_unlikely(new_lb != lb)) {
x                   LIST_REMOVE(p, entries);
x                   LIST_INSERT_HEAD(&moved_list_head, p, entries);
x               } 

                // indicate that this particle has been processed
                p->last_processed_time_ns = time_ns;
                __sync_fetch_and_add(&stats.processed,1);
                    
            }

            // unlock the particle_list_spinlock
            pthread_spin_unlock(&lb->particle_list_spinlock);

            // process the moved_list
            while ((p = moved_list_head.lh_first) != NULL) {
                LIST_REMOVE(p, entries);

                // determine the new locbox that the particle is being moved to
                new_lb = get_locbox(p->x_nm, p->y_nm, p->z_nm);

x              // if the new locbox is contained in a different radius shell then
x               // update the old and new radius shells
x               if (new_lb->radius_idx != lb->radius_idx) {
x                   radius_t * r_old = &radius[lb->radius_idx];
x                   if (!p->ion) {
x                       __sync_fetch_and_sub(&r_old->number_of_atoms,1);
x                   } else {
x                       __sync_fetch_and_sub(&r_old->number_of_ions,1);
x                   }
x                   __sync_fetch_and_sub(&r_old->sum_velocity_squared_m2pers2, p->velocity_squared_m2pers2);

x                   radius_t * r_new = &radius[new_lb->radius_idx];
x                   if (!p->ion) {
x                       __sync_fetch_and_add(&r_new->number_of_atoms,1);
x                   } else {
x                       __sync_fetch_and_add(&r_new->number_of_ions,1);
x                   }
x                   __sync_fetch_and_add(&r_new->sum_velocity_squared_m2pers2, p->velocity_squared_m2pers2);

x                   __sync_fetch_and_add(&stats.moved_radius,1);
x               }

                // add the particle to the particle_list in the new locbox
                pthread_spin_lock(&new_lb->particle_list_spinlock);
                LIST_INSERT_HEAD(&new_lb->particle_list_head, p, entries);
                pthread_spin_unlock(&new_lb->particle_list_spinlock);

                // stats update
                __sync_fetch_and_add(&stats.moved_locbox,1);
            }
        }

        // indicate completed to the control_thread
        pthread_barrier_wait(&barrier);
    }
}
#endif
