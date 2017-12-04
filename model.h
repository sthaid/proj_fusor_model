#ifndef __MODEL_H__
#define __MODEL_H__

#include <sys/queue.h>

//
// DEFINES
//

// params range
#define MIN_CHAMBER_RADIUS          0.050            // 1.97 inches
#define MAX_CHAMBER_RADIUS          0.100            // 3.93 inches
#define MIN_GRID_RADIUS             0.005            // 0.196 inches
#define MAX_GRID_RADIUS(chrad)      ((chrad) * .75) 
#define MIN_CHAMBER_PRESSURE        0.13332237       // 1 mtorr
#define MAX_CHAMBER_PRESSURE        133.32237        // 1000 mtorr
#define MAX_GRID_VOLTAGE            -50000.0
#define MAX_GRID_CURRENT            0.1

// cell
#define CELL_SIZE                   0.001
#define CELL_VOLUME                 1e-9
#define MAX_CELL                    200       // 2 * MAX_CHAMBER_RADIUS / CELL_SIZE

// shell
#define SHELL_SIZE                  0.001
#define MAX_SHELL                   100       // MAX_CHAMBER_RADIUS / SHELL_SIZE

// particles
#define AVG_SIM_PARTICLES_PER_CELL  10. // 10.0   XXX
#define MAX_CHAMBER_CELLS           4200000   // (4/3 * PI * MAX_CHAMBER_RADIUS^3) / CELL_VOLUME
#define MAX_PARTICLES               42000000  // MAX_CHAMBER_CELLS * AVG_SIM_PARTICLES_PER_CELL)

//
// TYPEDEFS
//

struct shell_s;

typedef struct {
    float chamber_radius;    // meters
    float grid_radius;       // meters
    float chamber_pressure;  // pascals
    float grid_voltage;      // volts
    float grid_current;      // amps
} params_t;

typedef struct particle_s {  // xxx organize
    LIST_ENTRY(particle_s) cell_entries;
    LIST_ENTRY(particle_s) work_entries;
    double time_last_processed;
    struct cell_s * cell;
    float x;
    float y;
    float z;
    float vx;
    float vy;
    float vz;
    float v;
    float v_squared;
    bool  ion;
} particle_t;

typedef struct cell_s {
    LIST_HEAD(head_s, particle_s) particle_list_head;
    struct shell_s * shell;
    float ax;
    float ay;
    float az;
} cell_t;

typedef struct shell_s {
    float   volume;
    int32_t number_of_atoms;
    int32_t number_of_ions;
    double  sum_v_squared;
    int64_t ionization_event_count;
    int64_t recombination_event_count;
    int64_t fusion_event_count;
    // XXX add list of cells
} shell_t;

//
// VARIABLES
//

params_t     params;

particle_t * particles;
int32_t      max_particles;

shell_t      shell[MAX_SHELL];
int32_t      max_shell;

cell_t       cell[MAX_CELL][MAX_CELL][MAX_CELL];

float        num_real_particles_in_cell;
float        num_real_particles_per_sim_particle;

double       time_model_secs;

//
// PROTOTYPES
//

void model_init(float chamber_radius, float grid_radius, float chamber_pressure,
        float grid_voltage, float grid_current);
void model_start(void);
void model_stop(void);
bool model_is_running(void);
void model_terminate(void);

#endif
