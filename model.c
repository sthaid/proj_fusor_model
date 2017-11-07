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
#include "util_misc.h"

// XXX todo
// - profiling
// - review all places fpu is used,such as sqrt etc.
// - shuffle the locbox_list to prevent lock contention

//
// defines
//

// model time increment
#define DELTA_T_NS  1L
#define DELTA_T (1e-9 * DELTA_T_NS)

// show the value of a define
#define SHOW_DEFINE(x) INFO("define %s = %s\n", #x, SHOW_DEFINE_STR(x))
#define SHOW_DEFINE_STR(x) #x

// for use in call to madvise
#define PAGE_SIZE 4096L
#define ROUND_UP(x,n) (((uint64_t)(x) + ((uint64_t)(n) - 1)) & ~((uint64_t)(n) - 1))

// locbox size in nanometers
#define LOCBOX_SIZE_NM  (LOCBOX_SIZE_MM * 1000000L)

// conversion functions  xxx need these in hdr file
#define METERS_TO_INCHES(m)     ((m) * 39.3701)
#define NM_TO_METERS(nm)        ((nm) * 1e-9)
#define METERS_TO_NM(m)         ((m) * 1e9)
#define NM_TO_MM(nm)            ((nm) * 1e-6)
#define MM_TO_NM(mm)            ((mm) * 1e6)
#define AMU_TO_KG(amu)          ((amu) * 1.66054e-27)
#define MTORR_TO_PASCAL(mtorr)  ((mtorr) * 0.13332237)
#define UTORR_TO_PASCAL(utorr)  ((utorr) * 0.00013332237)
#define F_TO_C(t)               (((t) - 32.0) * (5.0 / 9.0))
#define F_TO_K(t)               (F_TO_C(t) + 273.15)

// room temperature
#define ROOM_TEMPERATURE_KELVIN (F_TO_K(70.0))

// deuterium   // XXX or use D AMU
#define D2_AMU 4.03
#define D_AMU 2.000  // XXX
#define D2_KG  AMU_TO_KG(D2_AMU)

// kinetic temperature
// reference: http://hyperphysics.phy-astr.gsu.edu/hbase/Kinetic/kintem.html
#define k 1.38066e-23  // Boltzmann constant Joules/Kelvin
#define TEMPERATURE_TO_VELOCITY(t,m) (sqrt((t) * (3. * k / (m))))
#define VELOCITY_TO_TEMPERATURE(v,m) (((m) / (3. * k)) * (v) * (v))

// ideal gas law
// http://hyperphysics.phy-astr.gsu.edu/hbase/Kinetic/idegas.html
#define NUMBER_OF_MOLECULES(p,v,t)        ((p) * (v) / (k * (t)))
#define NUMBER_DENSITY_OF_MOLECULES(p,t)  ((p) / (k * (t)))

// Coulomb's Law
// http://hyperphysics.phy-astr.gsu.edu/hbase/electric/elefor.html
//    F = ke q1 * q2 / r^2
#define ke  8.987552e9  // Coulomb's Constant  N m^2 C^-2
#define COULOMB_FORCE(q1,q2,r)  ((ke * (q1) * (q2)) / ((r) * (r)))
#define ELECTRON_CHARGE (-1.60217662e-19)   // C
#define PROTON_CHARGE   (-ELECTRON_CHARGE)  // C

// http://physics.bu.edu/~duffy/PY106/Potential.html
//  V = ke Q / r
//  Q = V * r / ke
// XXX defines  ELECTRIC_POTENTIAL_TO_CHARGE ....


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

static int32_t           roomtemp_velocity_nmperdt;
static int32_t           roomtemp_velocity_mpers;
 
//
// prototypes
//

void init_particle(particle_t * p);
void * control_thread(void * cx);
void * work_thread(void * cx);

//
// inline functions
//

// XXX review
inline locbox_t * get_locbox(int32_t x_nm, int32_t y_nm, int32_t z_nm)
{
    int32_t x_idx = (x_nm + (MAX_LOCBOX*LOCBOX_SIZE_NM/2)) / LOCBOX_SIZE_NM;  // xxx use >>
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
    int64_t num_real_particles_in_locbox;
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
          METERS_TO_INCHES(params.chamber_radius_nm / 1e9));
    DEBUG("params.grid_radius_nm         = %d (%f inches)\n", 
          params.grid_radius_nm,
          METERS_TO_INCHES(params.grid_radius_nm / 1e9));
    DEBUG("params.chamber_pressure_utorr = %d\n", params.chamber_pressure_utorr);
    DEBUG("params.grid_voltage_v         = %d\n", params.grid_voltage_v);
    DEBUG("params.grid_current_ua        = %d\n", params.grid_current_ua);

    // initialize roomtemp_velocity_mpers and roomtemp_velocity_nmperdt
    roomtemp_velocity_mpers = TEMPERATURE_TO_VELOCITY(ROOM_TEMPERATURE_KELVIN, D2_KG);
    roomtemp_velocity_nmperdt = roomtemp_velocity_mpers * DELTA_T_NS;
    DEBUG("ROOM_TEMPERATURE_KELVIN    = %lf\n", ROOM_TEMPERATURE_KELVIN);
    DEBUG("roomtemp_velocity_mpers    = %d\n", roomtemp_velocity_mpers);
    DEBUG("roomtemp_velocity_nmperdt  = %d\n", roomtemp_velocity_nmperdt);

    // init max_radius
    max_radius = params.chamber_radius_nm / 1000000L / RADIUS_SHELL_SIZE_MM;
    DEBUG("max_radius      = %d\n", max_radius);

    // init locbox 
    q_grid = params.grid_voltage_v * NM_TO_METERS(params.grid_radius_nm) / ke;  // XXX macro
    for (x_idx = 0; x_idx < MAX_LOCBOX; x_idx++) {
        for (y_idx = 0; y_idx < MAX_LOCBOX; y_idx++) {
            for (z_idx = 0; z_idx < MAX_LOCBOX; z_idx++) {
                locbox_t * lb = &locbox[x_idx][y_idx][z_idx];

                LIST_INIT(&lb->particle_list_head);
                pthread_spin_init(&lb->particle_list_spinlock, PTHREAD_PROCESS_PRIVATE);

                x_nm = x_idx * LOCBOX_SIZE_NM - (MAX_LOCBOX*LOCBOX_SIZE_NM/2) + LOCBOX_SIZE_NM/2;
                y_nm = y_idx * LOCBOX_SIZE_NM - (MAX_LOCBOX*LOCBOX_SIZE_NM/2) + LOCBOX_SIZE_NM/2;
                z_nm = z_idx * LOCBOX_SIZE_NM - (MAX_LOCBOX*LOCBOX_SIZE_NM/2) + LOCBOX_SIZE_NM/2;
                lb->r_nm = hypotenuse(x_nm, y_nm, z_nm);
                lb->radius_idx = lb->r_nm / (RADIUS_SHELL_SIZE_MM * 1000000L);
                if (lb->radius_idx >= MAX_RADIUS) {
                    FATAL("radius_idx %d too big MAX_RADIUS = %ld\n", lb->radius_idx, MAX_RADIUS);
                }

                // XXX
                // COULOMB_FORCE(q1,q2,r)
                // q_grid = grid_voltage * grid_radius / ke
                // AMU_TO_KG(D2_AMU)

                if (lb->r_nm > params.grid_radius_nm) {
                    double f = COULOMB_FORCE(q_grid, PROTON_CHARGE, NM_TO_METERS(lb->r_nm));
                    double m = AMU_TO_KG(D_AMU);
                    lb->dxv_nmperdt = -(f / m) * ((double)x_nm / lb->r_nm) * DELTA_T * (1e9 * DELTA_T);
                    lb->dyv_nmperdt = -(f / m) * ((double)y_nm / lb->r_nm) * DELTA_T * (1e9 * DELTA_T);
                    lb->dzv_nmperdt = -(f / m) * ((double)z_nm / lb->r_nm) * DELTA_T * (1e9 * DELTA_T);
                } else {
                    lb->dxv_nmperdt = 0;
                    lb->dyv_nmperdt = 0;
                    lb->dzv_nmperdt = 0;
                }
                //if (y_idx == MAX_LOCBOX/2 && z_idx == MAX_LOCBOX/2) {
                    //INFO("x_nm = %d r_nm = %d  =>  %lf %lf %lf\n",
                            //x_nm, r_nm, dxv_mpersec, dxv_nmpersec, dxv_nmperdt);
                //}
            }
        }
    }

    // init locbox_list  
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

    // determine num_real_particles_per_virtual_particle 
    // using the ideal gas law
    num_real_particles_in_locbox = 
        NUMBER_OF_MOLECULES(UTORR_TO_PASCAL(params.chamber_pressure_utorr),
                            LOCBOX_VOLUME_CU_MM/1000000000.0,
                            ROOM_TEMPERATURE_KELVIN);
    num_real_particles_per_virtual_particle = num_real_particles_in_locbox / AVERAGE_PARTICLES_PER_LOCBOX;
    DEBUG("num_real_particles_in_locbox            = %ld\n", num_real_particles_in_locbox);
    DEBUG("num_real_particles_per_virtual_particle = %ld\n", num_real_particles_per_virtual_particle);

    // init time_ns
    time_ns = 0;

    // create threads
    max_worker_threads = sysconf(_SC_NPROCESSORS_ONLN) - 1;  // xxx comment save 1 for display
    DEBUG("creating control_thread and %d worker_threads ...\n", max_worker_threads);
    pthread_barrier_init(&barrier, NULL, max_worker_threads+1);
    pthread_create(&thread_id, NULL, control_thread, NULL);
    for (i = 0; i < max_worker_threads; i++) {
        pthread_create(&thread_id, NULL, work_thread, (void*)(intptr_t)i);
    }
    DEBUG("done creating threads\n");
}

void init_particle(particle_t * p)
{
    int32_t x_nm, y_nm, z_nm;
    int32_t x_dir, y_dir, z_dir, hypot, xv_nmperdt, yv_nmperdt, zv_nmperdt;
    locbox_t * lb;

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
    xv_nmperdt = (int64_t)x_dir * roomtemp_velocity_nmperdt / hypot;
    yv_nmperdt = (int64_t)y_dir * roomtemp_velocity_nmperdt / hypot;
    zv_nmperdt = (int64_t)z_dir * roomtemp_velocity_nmperdt / hypot;

#if 1
    // sanity check
    int32_t check_velocity_nmperdt = hypotenuse(xv_nmperdt, yv_nmperdt, zv_nmperdt);
    if (abs(check_velocity_nmperdt-roomtemp_velocity_nmperdt) > 2) {
        FATAL("roomtemp_velocity_nmperdt = %d  check_velocity_nmperdt = %d\n", 
              roomtemp_velocity_nmperdt, check_velocity_nmperdt);
    }
#endif

    // position and veloicty
    memset(p,0,sizeof(particle_t));
    p->x_nm = x_nm;
    p->y_nm = y_nm;
    p->z_nm = z_nm;
    p->xv_nmperdt = xv_nmperdt;
    p->yv_nmperdt = yv_nmperdt;
    p->zv_nmperdt = zv_nmperdt;
    p->velocity_squared_m2pers2 = (int64_t)roomtemp_velocity_mpers * (int64_t)roomtemp_velocity_mpers;

    // XXX
    p->ion = (lb->radius_idx == 63);

    // add to locbox list
    LIST_INSERT_HEAD(&lb->particle_list_head, p, entries);

    // updata associated radius info
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

// -----------------  MODEL THREADS  ---------------------------------------------------------

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
        // XXX try prefetch to improve performance
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

                // XXX velocity
                if (p->ion) {
                    p->xv_nmperdt += lb->dxv_nmperdt;
                    p->yv_nmperdt += lb->dyv_nmperdt;
                    p->zv_nmperdt += lb->dzv_nmperdt;
                }

                // update the particle's position, and 
                // determine the locbox it is now in following the position update
                p->x_nm += p->xv_nmperdt;
                p->y_nm += p->yv_nmperdt;
                p->z_nm += p->zv_nmperdt;
                new_lb = get_locbox(p->x_nm, p->y_nm, p->z_nm);

                // if the particle is now outside the chamber then 
                // - put the particle back where it was, and
                // - update its velocity to be roomtemp toward the center of the chamber
                if (__glibc_unlikely(new_lb->radius_idx >= max_radius)) {
                    // put the particle back where it was
                    p->x_nm -= p->xv_nmperdt;
                    p->y_nm -= p->yv_nmperdt;
                    p->z_nm -= p->zv_nmperdt;
                    new_lb = lb;

                    // update its velocity to be roomtemp toward the center of the chamber
                    p->xv_nmperdt = (int64_t)(-p->x_nm) * roomtemp_velocity_nmperdt / lb->r_nm;
                    p->yv_nmperdt = (int64_t)(-p->y_nm) * roomtemp_velocity_nmperdt / lb->r_nm;
                    p->zv_nmperdt = (int64_t)(-p->z_nm) * roomtemp_velocity_nmperdt / lb->r_nm;

#if 0 // xxx debug, delete
                    static bool first = true;
                    if (first) {
                        first = false;
                        int32_t vel = hypotenuse(p->xv_nmperdt, p->yv_nmperdt, p->zv_nmperdt);
                        DEBUG("vel = %d  roomtemp_vel = %d   vel_vect= %d %d %d  pos = %d %d %d\n",
                              vel, roomtemp_velocity_nmperdt, 
                              p->xv_nmperdt, p->yv_nmperdt, p->zv_nmperdt,
                              p->x_nm, p->y_nm, p->z_nm);
                    }
#endif
                }

                // if the particle has moved to another locbox then
                // remove it from this locbox list, and add it to the
                // moved_list; the reason for the moved_list is that we
                // can not put this particle on the new_lb particle_list 
                // at this time due to a possible deadlock; we first 
                // need to drop the particle_list_spinlock
                if (__glibc_unlikely(new_lb != lb)) {
                    LIST_REMOVE(p, entries);
                    LIST_INSERT_HEAD(&moved_list_head, p, entries);
                } 

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

                // if the new locbox is contained in a different radius shell then
                // update the old and new radius shells
                if (new_lb->radius_idx != lb->radius_idx) {
                    radius_t * r_old = &radius[lb->radius_idx];
                    if (!p->ion) {
                        __sync_fetch_and_sub(&r_old->number_of_atoms,1);
                    } else {
                        __sync_fetch_and_sub(&r_old->number_of_ions,1);
                    }
                    __sync_fetch_and_sub(&r_old->sum_velocity_squared_m2pers2, p->velocity_squared_m2pers2);

                    radius_t * r_new = &radius[new_lb->radius_idx];
                    if (!p->ion) {
                        __sync_fetch_and_add(&r_new->number_of_atoms,1);
                    } else {
                        __sync_fetch_and_add(&r_new->number_of_ions,1);
                    }
                    __sync_fetch_and_add(&r_new->sum_velocity_squared_m2pers2, p->velocity_squared_m2pers2);

                    __sync_fetch_and_add(&stats.moved_radius,1);
                }

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
