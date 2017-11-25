#ifndef __MODEL_H__
#define __MODEL_H__

// XXX 
// - should chamber pressure param be pascals

#include <sys/queue.h>

//
// DEFINES
//

// general
#define MB 0x100000L
#define CUBED(x) ((x) * (x) * (x))

// params range
#define MIN_CHAMBER_DIAMETER_MM     100    // 3.93 inches
#define MAX_CHAMBER_DIAMETER_MM     200    // 7.87 inches
#define MIN_GRID_DIAMETER_MM        10     // 0.39 inches
#define MAX_GRID_DIAMETER_MM(chdia) ((chdia) * 3 / 4) 
#define MIN_CHAMBER_PRESSURE_MTORR  1
#define MAX_CHAMBER_PRESSURE_MTORR  1000
#define MAX_GRID_VOLTAGE_KV         -100
#define MAX_GRID_CURRENT_MA         1000

// cell
#define CELL_SIZE_MM                1
#define CELL_SIZE_NM                (CELL_SIZE_MM * 1000000)
#define CELL_VOLUME_CU_MM           CUBED(CELL_SIZE_MM)
#define MAX_CELL                    (MAX_CHAMBER_DIAMETER_MM / CELL_SIZE_MM)

// shell
#define SHELL_SIZE_MM               1
#define SHELL_SIZE_NM               (SHELL_SIZE_MM * 1000000)
#define MAX_SHELL                   (MAX_CHAMBER_DIAMETER_MM / 2 / SHELL_SIZE_MM)

// particles
#define AVG_SIM_PARTICLES_PER_CELL  10
#define MAX_CHAMBER_VOLUME_CU_MM    ((int32_t)(CUBED(MAX_CHAMBER_DIAMETER_MM/2) * 4L * 3141593 / 3000000))
#define MAX_PARTICLES               (MAX_CHAMBER_VOLUME_CU_MM / CELL_VOLUME_CU_MM * AVG_SIM_PARTICLES_PER_CELL)

//
// TYPEDEFS
//

struct shell_s;

typedef struct {
    int32_t chamber_radius_nm;
    int32_t grid_radius_nm;
    int32_t chamber_pressure_utorr;
    int32_t grid_voltage_v;
    int32_t grid_current_ua;
} params_t;

typedef struct particle_s {
    LIST_ENTRY(particle_s) entries;
    int32_t x_nm;
    int32_t y_nm;
    int32_t z_nm;
    int32_t xv_nmperns;
    int32_t yv_nmperns;
    int32_t zv_nmperns;
    int64_t velocity_squared_m2pers2;
    bool    ion;
} particle_t;

typedef struct cell_s {
    LIST_HEAD(head_s, particle_s) particle_list_head;
    struct shell_s * shell;
    int32_t r_nm;
    int32_t xa_nmperns2;
    int32_t ya_nmperns2;
    int32_t za_nmperns2;
} cell_t;

typedef struct shell_s {
    int32_t volume_cu_mm;
    int32_t number_of_atoms;
    int32_t number_of_ions;
    int64_t ionization_event_count;
    int64_t recombination_event_count;
    int64_t fusion_event_count;
    int64_t sum_velocity_squared_m2pers2;
    // XXX add list of cells
} shell_t;

//
// VARIABLES
//

params_t     params;

particle_t   particles[MAX_PARTICLES];
int32_t      max_particles;

shell_t      shell[MAX_SHELL];
int32_t      max_shell;

cell_t       cell[MAX_CELL][MAX_CELL][MAX_CELL];

double       num_real_particles_per_sim_particle;

int64_t      time_ns;

//
// PROTOTYPES
//

void model_init_from_params(char * params_str);
void model_init_from_file(char * filename_str);
void model_start(void);
void model_stop(void);
void model_terminate(void);

#endif
