//
// TDOO
// - check the size of the particle and loc struct
//

//
// UNITS
//

// to optimize performance this model uses integer math, 
// these units are chosen with that in mind
//
// distance:  nanometers    
// velocity:  micrometers/sec  
// time       nanoseconds   
// pressure:  millipascals  
// voltage:   volts
// current:   microamps     

//
// DEFINES
//

#define MAX_PARTICLE          10000000L   // ten million
#define MAX_CHAMBER_DIAMETER  254000000L  // .254 m or 10 in
#define LOCATION_BOX_SIZE     1000000L    // 1 millimeter
#define RADIUS_SHELL_SIZE     1000000L    // 1 millimeter

#define MAX_LOC_XYZ (MAX_CHAMBER_DIAMETER/LOCATION_BOX_SIZE)
#define MAX_RADIUS  (MAX_CHAMBER_DIAMETER/2/RADIUS_SHELL_SIZE)


//
// TYPEDEFS
//

typedef struct {
    int64_t chamber_diameter;
    int64_t chamber_pressure;
    int64_t grid_diameter;
    int64_t grid_voltage;
    int64_t grid_current;
    int64_t delta_t;
} param_t;

typedef struct particle_s {
    int64_t x;
    int64_t y;
    int64_t z;
    int64_t r;
    int64_t x_velocity;
    int64_t y_velocity;
    int64_t z_velocity;
    int64_t velocity_squared;
    bool    ion;
    LIST_ENTRY(particle_s) entries;
} particle_t;

typedef struct {
    LIST_HEAD(head_s, particle_s) particle_list_head;
    bool inside_chamber;
} location_t;

typedef struct {
    int64_t temperature;
    int64_t atom_number_density;
    int64_t ion_number_density;
    int64_t ionization_count;
    int64_t recombination_count;
} radius_t;

//
// VARIABLES
//

param_t    param;

particle_t particle[MAX_PARTICLE];
location_t location[MAX_LOC_XYZ][MAX_LOC_XYZ][MAX_LOC_XYZ];
radius_t   radius[MAX_RADIUS];
int64_t    time;

