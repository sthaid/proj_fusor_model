#ifndef __MODEL_H__
#define __MODEL_H__

//
// XXX 
// - review use of int64 vs 32
// - check the size of the data allocations
// - CUBED eval once
// - should the defines really be longs
// - should chamber pressure param be pascals
// - num_real_particles_per_virtual_particle may be better as double
//

#include <sys/queue.h>

//
// DEFINES
//

// general
#define MB 0x100000L
#define CUBED(x) ((x) * (x) * (x))

// params range
#define MIN_CHAMBER_DIAMETER_MM       100L    // 3.93 inches
#define MAX_CHAMBER_DIAMETER_MM       200L    // 7.87 inches
#define MIN_GRID_DIAMETER_MM          10L     // 0.39 inches
#define MAX_GRID_DIAMETER_MM(chdia)   ((chdia) * 3L / 4L) 
#define MIN_CHAMBER_PRESSURE_MTORR    1L
#define MAX_CHAMBER_PRESSURE_MTORR    1000L
#define MAX_GRID_VOLTAGE_KV           100L
#define MAX_GRID_CURRENT_MA           1000L

// locbox
#define LOCBOX_SIZE_MM                1L
#define LOCBOX_SIZE_NM                (LOCBOX_SIZE_MM * 1000000L)
#define LOCBOX_VOLUME_CU_MM           CUBED(LOCBOX_SIZE_MM)
#define MAX_LOCBOX                    ((MAX_CHAMBER_DIAMETER_MM + 10L) / LOCBOX_SIZE_MM)
 // XXX ^^^ should this be even number

// radius
#define RADIUS_SHELL_SIZE_MM          1L   // must be >= LOCBOX_SIZE_MM
#define MAX_RADIUS                    (300L / RADIUS_SHELL_SIZE_MM)

// particles
#define AVERAGE_PARTICLES_PER_LOCBOX  1L  // XXX was 10L
#define MAX_CHAMBER_VOLUME_CU_MM      (CUBED(MAX_CHAMBER_DIAMETER_MM/2) * 4L * 3141593L / 3000000L) 
#define MAX_PARTICLES                 (MAX_CHAMBER_VOLUME_CU_MM / LOCBOX_VOLUME_CU_MM * AVERAGE_PARTICLES_PER_LOCBOX)

//
// TYPEDEFS
//

typedef struct {
    int32_t chamber_radius_nm;
    int32_t grid_radius_nm;
    int32_t chamber_pressure_utorr;
    int32_t grid_voltage_v;
    int32_t grid_current_ua;
} params_t;

typedef struct particle_s {
    LIST_ENTRY(particle_s) entries;
    int64_t last_processed_time_ns;
    int32_t x_nm;
    int32_t y_nm;
    int32_t z_nm;
    int32_t xv_nmperdt;
    int32_t yv_nmperdt;
    int32_t zv_nmperdt;
    int64_t velocity_squared_m2pers2;
    bool    ion;
} particle_t;

typedef struct locbox_s {
    LIST_HEAD(head_s, particle_s) particle_list_head;
    pthread_spinlock_t particle_list_spinlock;
    int32_t radius_idx;
    int32_t r_nm;
    int32_t dxv_nmperdt;
    int32_t dyv_nmperdt;
    int32_t dzv_nmperdt;
} locbox_t;

typedef struct {
    int64_t volume_cu_mm;
    int64_t number_of_atoms;
    int64_t number_of_ions;
    int64_t ionization_event_count;
    int64_t recombination_event_count;
    int64_t fusion_event_count;
    int64_t sum_velocity_squared_m2pers2;
} radius_t;

//
// VARIABLES
//

params_t     params;

particle_t   particles[MAX_PARTICLES];
int32_t      max_particles;

radius_t     radius[MAX_RADIUS];
int32_t      max_radius;

locbox_t     locbox[MAX_LOCBOX][MAX_LOCBOX][MAX_LOCBOX];

int64_t      time_ns;
int64_t      num_real_particles_per_virtual_particle;

//
// PROTOTYPES
//

void model_init_from_params(char * params_str);
void model_init_from_file(char * filename_str);
void model_start(void);
void model_stop(void);
void model_terminate(void);

#endif
