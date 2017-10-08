#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include <assert.h>

#include <pthread.h>
#include <sys/mman.h>

#include "model.h"
#include "util_misc.h"

// XXX todo
// - profiling
// - review all places fpu is used,such as sqrt etc.

//
// defines
//

// XXX review this section
#define LOCBOX_SIZE_NM                (1000000L * LOCBOX_SIZE_MM)
#define DELTA_T_NS  1L

#define PAGE_SIZE 4096L
#define ROUND_UP(x,n) (((uint64_t)(x) + ((uint64_t)(n) - 1)) & ~((uint64_t)(n) - 1))

#define INCHES_PER_METER  39.3701

#define D2_AMU 4.03                // Deuterium molecule mass
#define D2_KG  AMU_TO_KG(D2_AMU)
#define AMU_TO_KG(amu)          ((amu) * 1.66054e-27)

// kinetic temperature
// reference: http://hyperphysics.phy-astr.gsu.edu/hbase/Kinetic/kintem.html
#define k 1.38066e-23  // Boltzmann constant Joules/Kelvin
#define TEMPERATURE_TO_VELOCITY(t,m) (sqrt((t) * (3. * k / (m))))
#define VELOCITY_TO_TEMPERATURE(v,m) (((m) / (3. * k)) * (v) * (v))


// 
// typedefs
//

//
// variables
//

static locbox_t        * locbox_list[sizeof(locbox)/sizeof(locbox_t)];
static int32_t           max_locbox_list;
static int32_t           locbox_list_idx;

static bool              run_request;
static bool              running;

static pthread_barrier_t barrier;
 
//
// prototypes
//

void init_particle(particle_t * p);
void * control_thread(void * cx);
void * work_thread(void * cx);

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

inline int32_t random_range(int32_t min, int32_t max)
{
    int64_t extent = (int64_t)max - min + 1L;
    return random() * extent / (RAND_MAX+1L) + min;
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
    int32_t x_nm, y_nm, z_nm, r_nm;
    int32_t max_worker_threads;
    pthread_t thread_id;

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
    if (grid_voltage_kv < 0 || 
        grid_voltage_kv > MAX_GRID_VOLTAGE_KV)  {
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
          params.chamber_radius_nm / 1e9 * INCHES_PER_METER);
    DEBUG("params.grid_radius_nm         = %d (%f inches)\n", 
          params.grid_radius_nm,
          params.grid_radius_nm / 1e9 * INCHES_PER_METER);
    DEBUG("params.chamber_pressure_utorr = %d\n", params.chamber_pressure_utorr);
    DEBUG("params.grid_voltage_v         = %d\n", params.grid_voltage_v);
    DEBUG("params.grid_current_ua        = %d\n", params.grid_current_ua);

    // init max_radius
    max_radius = params.chamber_radius_nm / 1000000L / RADIUS_SHELL_SIZE_MM;
    DEBUG("max_radius      = %d\n", max_radius);

    // init locbox 
    for (x_idx = 0; x_idx < MAX_LOCBOX; x_idx++) {
        for (y_idx = 0; y_idx < MAX_LOCBOX; y_idx++) {
            for (z_idx = 0; z_idx < MAX_LOCBOX; z_idx++) {
                locbox_t * lb = &locbox[x_idx][y_idx][z_idx];

                LIST_INIT(&lb->particle_list_head);
                pthread_spin_init(&lb->particle_list_spinlock, PTHREAD_PROCESS_PRIVATE);

                x_nm = x_idx * LOCBOX_SIZE_NM - (MAX_LOCBOX*LOCBOX_SIZE_NM/2) + LOCBOX_SIZE_NM/2;
                y_nm = y_idx * LOCBOX_SIZE_NM - (MAX_LOCBOX*LOCBOX_SIZE_NM/2) + LOCBOX_SIZE_NM/2;
                z_nm = z_idx * LOCBOX_SIZE_NM - (MAX_LOCBOX*LOCBOX_SIZE_NM/2) + LOCBOX_SIZE_NM/2;
                r_nm = hypotenuse(x_nm, y_nm, z_nm);
                lb->radius_idx = r_nm / (RADIUS_SHELL_SIZE_MM * 1000000L);
                if (lb->radius_idx >= MAX_RADIUS) {
                    FATAL("radius_idx %d too big MAX_RADIUS = %ld\n", lb->radius_idx, MAX_RADIUS);
                }
            }
        }
    }

    // init locbox_list  XXX shuffle 
    for (x_idx = 0; x_idx < MAX_LOCBOX; x_idx++) {
        for (y_idx = 0; y_idx < MAX_LOCBOX; y_idx++) {
            for (z_idx = 0; z_idx < MAX_LOCBOX; z_idx++) {
                locbox_t * lb = &locbox[x_idx][y_idx][z_idx];
                if (lb->radius_idx < max_radius) {
                    locbox_list[max_locbox_list++] = &locbox[x_idx][y_idx][z_idx];
                }
            }
        }
    }
    DEBUG("max_locbox_list = %d\n", max_locbox_list);

    // init radius
    for (x_idx = 0; x_idx < MAX_LOCBOX; x_idx++) {
        for (y_idx = 0; y_idx < MAX_LOCBOX; y_idx++) {
            for (z_idx = 0; z_idx < MAX_LOCBOX; z_idx++) {
                locbox_t * lb = &locbox[x_idx][y_idx][z_idx];
                radius[lb->radius_idx].volume_cu_mm += LOCBOX_VOLUME_CU_MM;
            }
        }
    }

    // init particles and max_particles
    max_particles = max_locbox_list * AVERAGE_PARTICLES_PER_LOCBOX;
    DEBUG("max_particles   = %d\n", max_particles);
    DEBUG("initializing particles ...\n");
    for (i = 0; i < max_particles; i++) {
        init_particle(&particles[i]);
    }
    DEBUG("done initializing particles\n");

    // create threads
    max_worker_threads = sysconf(_SC_NPROCESSORS_ONLN);
    DEBUG("creating control_thread and %d worker_threads ...\n", max_worker_threads);
    pthread_barrier_init(&barrier, NULL, max_worker_threads+1);
    pthread_create(&thread_id, NULL, control_thread, NULL);
    for (i = 0; i < max_worker_threads; i++) {
        pthread_create(&thread_id, NULL, work_thread, (void*)(intptr_t)i);
    }
    DEBUG("done creating threads\n");

    // init time_ns
    time_ns = 0;

    // XXX
    // determine num_real_particles_per_virtual_particle based
    // on the pressure param
    // 
    // num_real_particles_in_locbox = P V / (R T)
    // num_real_particles_per_virtual_particle = num_real_particles_in_locbox / 
    //                                           VIRT_PARTICLES_IN_LOCBOX
    //   where
    //     P is pressure param
    //     V is volume of a locbox element
    //     T is room temperature
}

void init_particle(particle_t * p)
{
    int32_t x_nm, y_nm, z_nm;
    int32_t x_dir, y_dir, z_dir, hypot, xv_nmperdt, yv_nmperdt, zv_nmperdt;
    locbox_t * lb;

    static int64_t velocity_nmperdt, velocity_mpers;

    // initialize velocity_nmperdt, if not already initialized;
    // all particles will be given this velocity which is equivalent to 300 degrees kelvin
    if (velocity_nmperdt == 0) {
        velocity_mpers = TEMPERATURE_TO_VELOCITY(300,AMU_TO_KG(D2_AMU));
        velocity_nmperdt = velocity_mpers * DELTA_T_NS;
        DEBUG("velocity_mpers    = %ld\n", velocity_mpers);
        DEBUG("velocity_nmperdt  = %ld\n", velocity_nmperdt);
    }

    // get a random location 
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
    // - the direction is random, the magnitude is 300 kelvin
    // - the random direction is determined by picking a random point within a sphere
    while (true) {
        x_dir = random_range(-1000000, 1000000);
        y_dir = random_range(-1000000, 1000000);
        z_dir = random_range(-1000000, 1000000);
        hypot = hypotenuse(x_dir,y_dir,z_dir);
        if (hypot < 1000000) {
            break;
        }
    }
    xv_nmperdt = (int64_t)x_dir * velocity_nmperdt / hypot;
    yv_nmperdt = (int64_t)y_dir * velocity_nmperdt / hypot;
    zv_nmperdt = (int64_t)z_dir * velocity_nmperdt / hypot;

#if 1
    // sanity check
    int64_t check_velocity_nmperdt = hypotenuse(xv_nmperdt, yv_nmperdt, zv_nmperdt);
    if (abs(check_velocity_nmperdt-velocity_nmperdt) > 2) {
        FATAL("velocity_nmperdt = %ld  check_velocity_nmperdt = %ld\n", 
              velocity_nmperdt, check_velocity_nmperdt);
    }
#endif

    // init particle: position and veloicty
    memset(p,0,sizeof(particle_t));
    p->x_nm = x_nm;
    p->y_nm = y_nm;
    p->z_nm = z_nm;
    p->xv_nmperdt = xv_nmperdt;
    p->yv_nmperdt = yv_nmperdt;
    p->zv_nmperdt = zv_nmperdt;
    p->velocity_squared_m2pers2 = velocity_mpers * velocity_mpers;

    // init particle: add to locbox list
    p->lb = lb;
    LIST_INSERT_HEAD(&lb->particle_list_head, p, entries);

    // init particle: updata associated radius info
    radius_t * r = &radius[lb->radius_idx];
    r->number_of_atoms++;
    r->sum_velocity_squared_m2pers2 +=  p->velocity_squared_m2pers2;
}

// -----------------  MODEL INIT FROM FILE  ----------------------------------------------

void model_init_from_file(char * filename_str)
{
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

// -----------------  MODEL ------------------------------------------------------------------

static uint64_t processed, moved;
void * control_thread(void * cx)
{
#if 0 // XXX
    uint64_t start_us, done_us;

    while (true) {
        // wait for run_flag to be set
        if (!run_request) {
            running = false;
            while (run_request == false) {
                usleep(1000);
            }
            running = true;
        }

        // advance time
        time_ns += DELTA_T_NS;

        // xxxx
        locbox_list_idx = 0;

        processed = 0;
        moved = 0;

        // release threads
        INFO("RELEASE THREADS for time %ld ns\n", time_ns);
        start_us = microsec_timer();
        pthread_barrier_wait(&barrier);

        // wait for threads to complete
        pthread_barrier_wait(&barrier);
        done_us = microsec_timer();
        INFO("THREADS DONE for time %ld ns - duration %ld us - processed %ld moved %ld\n", 
             time_ns, done_us - start_us,
             processed, moved);
        //sleep(5); //XXX
    }
#endif
    return NULL;
}

void * work_thread(void * cx)
{
#if 0 // XXX
    int32_t    thread_id  __attribute__ ((unused)) = (int64_t)cx;
    int32_t    idx;
    locbox_t * lb;
    particle_t * p;

    LIST_HEAD(head_s, particle_s) moved_list_head;

    while (true) {
        // wait to be started 
        pthread_barrier_wait(&barrier);
        //INFO("starting thread_id %d\n", thread_id);

        // loop until no more location boxes remain
        // - get the next location box 
        // - process all the particles within the location box
        //   that have not already been processed 
        while (true) {
            // get the next location box, if no more then break
            idx = __sync_fetch_and_add(&locbox_list_idx, 1);
            if (idx >= max_locbox_list) {
                break;
            }
            lb = locbox_list[idx];

            // AAA locking
            pthread_spin_lock(&lb->particle_list_spinlock);
    
            LIST_INIT(&moved_list_head);

            // loop over all particles within the location box,
            // and process them if they have not already been
            // processed; a particle may already have been processed
            // if it has just been moved into this location box
            particle_t * next_p;
            for (p = lb->particle_list_head.lh_first; p != NULL; p = next_p) {
                next_p = p->entries.le_next;
                //locbox_t * xxxlb, * afterlb;
                //xxxlb = get_locbox(p->x_nm, p->y_nm, p->z_nm);
                //INFO("lb %ld   xxxlb %ld\n",
                    //lb - &locbox[0][0][0],
                    //xxxlb - &locbox[0][0][0]);

                if (p->last_processed_time_ns == time_ns) {
                    continue;
                }

                // AAA acceleration of ions, store the accel in locbox

                //INFO("p_idx  %ld\n", p - particles);
                //INFO("before %d %d %d\n", p->x_nm, p->y_nm, p->z_nm);
                p->x_nm += p->xv_nmperdt;
                p->y_nm += p->yv_nmperdt;
                p->z_nm += p->zv_nmperdt;
                //INFO("after %d %d %d\n", p->x_nm, p->y_nm, p->z_nm);
                //afterlb = get_locbox(p->x_nm, p->y_nm, p->z_nm);
                //INFO("lb %ld   afterlb %ld\n",
                    //lb - &locbox[0][0][0],
                    //afterlb - &locbox[0][0][0]);

                locbox_t * new_lb;
                new_lb = get_locbox(p->x_nm, p->y_nm, p->z_nm);

                // XXX if outside chamber then put it back , and set new dir
                if (new_lb->radius_idx >= max_radius) {
                    //INFO("outside chamber ridx %d\n", new_lb->radius_idx);
                    p->x_nm -= p->xv_nmperdt;
                    p->y_nm -= p->yv_nmperdt;
                    p->z_nm -= p->zv_nmperdt;
                    new_lb = lb;

                    // AAA set new dir towards center of chamber
                    int64_t r_nm = hypotenuse(p->x_nm, p->y_nm, p->z_nm);  // AAA use lb->xxxx
                    int64_t target_velocity_mpers = 1000;  // XXX calc from physics
                    int64_t target_velocity_nmperdt = target_velocity_mpers * DELTA_T_NS;

                    p->xv_nmperdt = -p->x_nm * target_velocity_nmperdt / r_nm;
                    p->yv_nmperdt = -p->y_nm * target_velocity_nmperdt / r_nm;
                    p->zv_nmperdt = -p->z_nm * target_velocity_nmperdt / r_nm;

                    static bool first = true;
                    if (first) {
                        first = false;
                        int32_t vel = hypotenuse(p->xv_nmperdt, p->yv_nmperdt, p->zv_nmperdt);
        
                        INFO("vel = %d  target_vel = %ld   vel_vect= %d %d %d  pos = %d %d %d\n",
                            vel, target_velocity_nmperdt, 
                            p->xv_nmperdt, p->yv_nmperdt, p->zv_nmperdt,
                            p->x_nm, p->y_nm, p->z_nm);
                    }
                }

                p->last_processed_time_ns = time_ns;

                // XXX needs locking
                __sync_fetch_and_add(&processed,1);
                if (new_lb != lb) {
                    // AAA
                    //INFO("moving to new locbox\n");
                    // remove from list
                    LIST_REMOVE(p, entries);
                    LIST_INSERT_HEAD(&moved_list_head, p, entries);
                    // add to other list
                    // XXX also need to flag so it is not reprocessed;
                    //     maybe set thetime idx of the particle
                } 

                    
            }

            // AAA locking
            pthread_spin_unlock(&lb->particle_list_spinlock);

            // AAA moved list
            radius_t * r;
            while ((p = moved_list_head.lh_first) != NULL) {
                LIST_REMOVE(moved_list_head.lh_first, entries);
                r = &radius[p->lb->radius_idx];
                __sync_fetch_and_sub(&r->number_of_atoms,1); // AAA ions
                __sync_fetch_and_sub(&r->sum_velocity_squared_m2pers2, p->velocity_squared_m2pers2);

                locbox_t * new_lb;
                new_lb = get_locbox(p->x_nm, p->y_nm, p->z_nm);
                pthread_spin_lock(&new_lb->particle_list_spinlock);
                LIST_INSERT_HEAD(&new_lb->particle_list_head, p, entries);
                p->lb = new_lb;
                pthread_spin_unlock(&new_lb->particle_list_spinlock);
                r = &radius[p->lb->radius_idx];
                __sync_fetch_and_add(&r->number_of_atoms,1); // AAA ions
                __sync_fetch_and_add(&r->sum_velocity_squared_m2pers2, p->velocity_squared_m2pers2);

                __sync_fetch_and_add(&moved,1);
            }
            // AAA radius stuff,  use atomics
        }

        // indicate completed
        //INFO("done thread_id %d\n", thread_id);
        pthread_barrier_wait(&barrier);
    }
#endif
    return NULL;
}
