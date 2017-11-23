// XXX try running w/o secondary and compare percent
//            with secondary 2.077
//            w/o secondary  2.02   velocity = v1
//            w/o secondary  2.05   velocity = v2
// XXX try running at higher voltage, and verify the initial step size
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

#include "physics.h"
#include "util_misc.h"

// to build for UNIT_TEST:
// gcc -g -O2 -Wall -Iutil -DUNIT_TEST -lm -o ut model_ionization.c util/util_misc.c

//
// defines
//

#define MAX_IONIZATION_COUNT 10000000

//
// typedefs
//

//
// variables
//

static double  g_chamber_radius;
static double  g_grid_radius;
static double  g_chamber_pressure;
static double  g_grid_voltage;
static double  g_grid_current;
static double  g_ionization_count_delta;
static int32_t g_ionization_count[MAX_IONIZATION_COUNT];

//
// prototypes
//

static void simulate_electron(double position, double velocity);
static double ionization_cross_section(double ev);

// --------------------------------------------------------------------------------------

void model_ionization(
    double chamber_radius,
    double grid_radius,
    double chamber_pressure,  // pascals
    double grid_voltage,
    double grid_current,
    double ionization_count_delta)
{
    double time_interval_secs, t;
    int32_t i;

    #define MAX 100000
    #define N   100000

    INFO("starting\n");

    // XXX check if results are cached in file .model_ionization.dat

    // save args in global variables
    g_chamber_radius         = chamber_radius;
    g_grid_radius            = grid_radius;
    g_chamber_pressure       = chamber_pressure;
    g_grid_voltage           = grid_voltage;
    g_grid_current           = grid_current;
    g_ionization_count_delta = ionization_count_delta;

    // electrons will be created at the grid, one at a time,
    // at a time interval equivalent to the grid current

    time_interval_secs = -ELECTRON_CHARGE / g_grid_current;
    INFO("time_interval = %e\n", time_interval_secs);

    t = 0;
    for (i = 0; i < MAX; i++) {
        simulate_electron(g_grid_radius, 0);
        t += time_interval_secs;

        if ((i % N) == N - 1) {
            int32_t sum = 0, j;
            for (j = 0; j < 80; j++) {
                INFO("RESULT  %3d %4d\n", j, g_ionization_count[j]);
                sum += g_ionization_count[j];
            }
            INFO("RESULT i=%d sum=%d percent=%f t=%f ns\n", 
                 i+1, sum, (double)sum/(i+1)*100, t*1e9);
            BLANK_LINE;
        }
    }

    for (i = 0; i < 80; i++) {
        INFO("XXX RESULT i=%d  %e\n",
            i, g_ionization_count[i] / t);
    }

    // XXXcalculate results, which are an array, indexed by radius, 
    // XXXof ionization rate in units count/sec

    // XXX print results

    // XXX save results to file .model_ionization.dat

}

static void simulate_electron(
    double position,  
    double velocity)
{
    double q, number_density, force, acceleration, ev, cross_section, mfp;
    double delta_position, delta_t, probability_of_zero_collision_events;
    bool   ionization_occurred;

    // init constants
    q = POTENTIAL_TO_CHARGE(g_grid_voltage,g_grid_radius);
    number_density = NUMBER_DENSITY(g_chamber_pressure,ROOM_TEMPERATURE_KELVIN);

    // this loop will execute until position >= g_chamber_radius
    while (true) {
        // determine the electron's acceleration
        force = COULOMB_FORCE(ELECTRON_CHARGE,q,position);
        acceleration = force / ELECTRON_MASS;    
        //INFO("force=%f  accel=%f\n", force, acceleration);

        // determine the ionization cross section
        ev = JOULES_TO_EV(.5*ELECTRON_MASS*velocity*velocity);
        cross_section = ionization_cross_section(ev);

        // if the cross_section is zero then the mean-free-path will be infinite; 
        // so there will be no more ionizations caused by this electron
        // XXX what about higher voltages
        // XXX fix the comment
        if (cross_section == 0) {
            //INFO("x=%f v=%f ev=%f cs=%f\n", 
             //    position, velocity, ev, cross_section);
            if (ev <= 15) {
                delta_t = 1e-12;  // XXX how to verify this is small enough
                position += velocity * delta_t;
                velocity += acceleration * delta_t;
                continue;
            } else {
                //INFO("XXXXXXX  RETURN  x=%f  ev=%f\n", position, ev);
                return;
            }
        }

        // determine the mean free path, this is the average distance the
        // electron travels between ionization events
        mfp = MEAN_FREE_PATH(cross_section,number_density);
        //INFO("position=%f  ev=%f  mfp=%f\n", position, ev, mfp);

        // set delta_t to the time to travel 1/10 mfp  XXX comment
        //delta_position = 0.1 * mfp;
        delta_position = 0.00001 * mfp;
        delta_t = delta_position / velocity;
        //INFO("x=%f v=%f ev=%f cs=%f mfp=%f dx=%f dt=%f\n", 
             //position, velocity, ev, cross_section, mfp, delta_position, delta_t);

        // using Poison Distribution determine the probability of zero 
        // collision events occurring over the delta_t interval
        //  
        // https://en.wikipedia.org/wiki/Poisson_distribution
        //
        //                                       lambda^k
        // P(k-events-in-interval) = e^-lambda * --------
        //                                         k!
        //
        // where: lamda = average number of events per interval
        //        k     = 0,1,2...
        //
        // lambda = 0.00001, k = 0  ==> P = .99999
        probability_of_zero_collision_events = .99999;

        // choose a random number and use it along with the calculated probability 
        // to choose whether or not an ionization event occurred
        ionization_occurred = ((double)random()/RAND_MAX > probability_of_zero_collision_events);

        // if ionization did not occure then
        //   adjust the electron's speed and velocity;
        //   if the electron has reached the chamber radius then return
        // else
        //   increment ionization count at the electron's current position
        //   create a new electron: determine the velocity of the new
        //    electron and the existing electron 
        //   simulate the new electron, and continue simulating the existing electron
        // endif
        if (ionization_occurred == false) {
            position = position + delta_position;
            velocity = velocity + acceleration * delta_t;
            if (position >= g_chamber_radius) {
                //INFO("=======  RETURN  x=%f  ev=%f\n", position, ev);
                return;
            }
        } else {
            int32_t idx;
            double ke;
            //INFO("XXX IONIZATION OCCURRED position=%f velocity=%f\n",
                 //position, velocity);

            // increment ionization count for the current position
            idx = position / g_ionization_count_delta;
            if (idx < MAX_IONIZATION_COUNT) {
                INFO("XXX idx %d\n", idx);
                g_ionization_count[idx]++;
            } else {
                ERROR("idx %d >= MAX_IONIZATION_COUNT\n", idx);
            }

            ke = KINETIC_ENERGY(ELECTRON_MASS,velocity);
            ke = ke - EV_TO_JOULES(HYDROGEN_IONIZATION_ENERGY_EV);
            velocity = KINETIC_ENERGY_TO_VELOCITY(ke,ELECTRON_MASS);

            simulate_electron(position, velocity);
            velocity = 0;
        }
    }
}


// returns cross section in m^2
static double ionization_cross_section(double ev) 
{
    // Reference:
    //
    // https://www.nist.gov/pml/electron-impact-cross-sections-ionization-and-excitation-database
    //    Table of Atoms -> H -> Total
    //
    // https://physics.nist.gov/cgi-bin/Ionization/ion_data.php?id=HI&ision=I&initial=&total=Y

    static struct {
        double ev;
        double cs;   // cross_section  1e-16 cm^2
    } tbl[] = {
        { 14.60,  0.04430  },
        { 14.80,  0.05350  },
        { 15.00,  0.06260  },
        { 15.10,  0.06720  },
        { 15.20,  0.07180  },
        { 15.40,  0.08090  },
        { 15.60,  0.09000  },
        { 15.90,  0.10360  },
        { 16.10,  0.11260  },
        { 16.40,  0.12600  },
        { 16.60,  0.13490  },
        { 16.90,  0.14800  },
        { 17.10,  0.15670  },
        { 17.40,  0.16950  },
        { 17.60,  0.17790  },
        { 17.90,  0.19040  },
        { 18.10,  0.19850  },
        { 18.40,  0.21060  },
        { 18.70,  0.22250  },
        { 19.00,  0.23410  },
        { 19.30,  0.24540  },
        { 19.60,  0.25650  },
        { 20.00,  0.27100  },
        { 20.40,  0.28500  },
        { 20.90,  0.30190  },
        { 21.40,  0.31820  },
        { 22.00,  0.33690  },
        { 22.60,  0.35470  },
        { 23.30,  0.37440  },
        { 24.00,  0.39300  },
        { 24.80,  0.41290  },
        { 25.60,  0.43150  },
        { 26.60,  0.45300  },
        { 27.30,  0.46700  },
        { 28.30,  0.48560  },
        { 29.30,  0.50260  },
        { 30.50,  0.52120  },
        { 31.60,  0.53650  },
        { 32.80,  0.55170  },
        { 34.10,  0.56630  },
        { 35.40,  0.57940  },
        { 36.70,  0.59090  },
        { 38.10,  0.60200  },
        { 39.60,  0.61220  },
        { 41.20,  0.62170  },
        { 42.90,  0.63010  },
        { 44.70,  0.63760  },
        { 46.60,  0.64400  },
        { 48.60,  0.64940  },
        { 50.70,  0.65360  },
        { 52.90,  0.65680  },
        { 55.20,  0.65900  },
        { 57.60,  0.66010  },
        { 60.10,  0.66030  },
        { 63.00,  0.65940  },
        { 66.00,  0.65750  },
        { 69.00,  0.65470  },
        { 72.10,  0.65110  },
        { 75.50,  0.64660  },
        { 79.50,  0.64050  },
        { 84.00,  0.63310  },
        { 89.00,  0.62420  },
        { 94.00,  0.61500  },
        { 102.00, 0.59980  },
        { 103.00, 0.59790  },
        { 113.00, 0.57890  },
        { 121.00, 0.56390  },
        { 130.20, 0.54710  },
        { 138.20, 0.53310  },
        { 148.20, 0.51620  },
        { 158.20, 0.50030  },
        { 168.20, 0.48510  },
        { 178.20, 0.47080  },
        { 188.20, 0.45720  },
        { 198.20, 0.44440  },
        { 213.20, 0.42650  },
        { 228.20, 0.40990  },
        { 248.20, 0.38970  },
        { 268.20, 0.37150  },
        { 288.00, 0.35510  },
        { 317.90, 0.33300  },
        { 347.90, 0.31360  },
        { 387.90, 0.29120  },
        { 427.90, 0.27190  },
        { 467.90, 0.25510  },
        { 508.20, 0.24030  },
        { 548.20, 0.22740  },
        { 598.20, 0.21310  },
        { 668.20, 0.19610  },
        { 700.00, 0.18930  },
        { 748.20, 0.17990  },
        { 818.20, 0.16790  },
        { 898.20, 0.15610  },
        { 998.20, 0.14360  },
        { 1100.00, 0.13290 },
        { 1200.00, 0.12400 },
        { 1300.00, 0.11620 },
        { 1506.70, 0.10310 },
        { 1662.70, 0.09510 },
        { 1848.10, 0.08710 },
        { 1998.10, 0.08170 },
        { 2198.10, 0.07540 },
        { 2448.10, 0.06890 },
        { 2698.10, 0.06350 },
        { 2998.10, 0.05810 },
        { 3298.10, 0.05350 },
        { 3648.10, 0.04910 },
        { 3998.10, 0.04540 },
        {50000.00, 0.0     }  // note: this entry not from nist.gov, I added this entry
            };                // to extrapolate beyond 3998 ev.

    #define MAX_TBL (sizeof(tbl)/sizeof(tbl[0]))
    #define MAX_EV 50000
    static double ev_to_cross_section[MAX_EV];
    static bool ev_to_cross_section_initialized;

    // convert nist.gov table to ev_to_cross_section which is indexed by ev
    if (!ev_to_cross_section_initialized) {
        int32_t ev, i;
        for (ev = 0; ev < MAX_EV; ev++) {
            // if ev is not in range covered by table then 
            // set new_ebl cross section entry to 0
            if (ev < tbl[0].ev || ev > tbl[MAX_TBL-1].ev) {
                ev_to_cross_section[ev] = 0;
                continue;
            }

            // find location in the nist.gov tbl to interpolate
            for (i = 0; i < MAX_TBL-1; i++) {
                if (ev >= tbl[i].ev && ev <= tbl[i+1].ev) {
                    break;
                }
            }
            assert(i < MAX_TBL-1);

            // determine cross section by interpolating from tbl  XXX check this
            double slope = (tbl[i+1].cs - tbl[i].cs) / (tbl[i+1].ev - tbl[i].ev);
            double cs = tbl[i].cs + slope * (ev - tbl[i].ev);    // 1e-16 cm^2

            // store cross section in ev_to_cross_section, in units m^2
            ev_to_cross_section[ev] = cs * 1e-20;  // m^2
        }

        // set initialized flag
        ev_to_cross_section_initialized = true;
    }

    // return cross section, in m^2 via table lookup
    if (ev >= MAX_EV) {
        return 0;
    }
    return ev_to_cross_section[(int32_t)ev];
}

// --------------------------------------------------------------------------------------

#ifdef UNIT_TEST

int32_t main(int32_t argc, char **argv)
{
    //double ke = KINETIC_ENERGY(ELECTRON_MASS, 100);
    //double v = KINETIC_ENERGY_TO_VELOCITY(ke,ELECTRON_MASS);
    //INFO("XXX v %f\n", v);
    //exit(1);

    model_ionization(
        //MM_TO_M(75.),       // chamber_radius, meters
        MM_TO_M(150.),       // chamber_radius, meters
        MM_TO_M(19.),       // grid_radius, meters
        MTORR_TO_PASCAL(2*15),  // chamber_pressure, pascals  XXX why 2
        -30000,                 // grid_voltage, volts
        .005,                   // grid_current, amps
        .001);                  // ionization_count_delta, meters


    INFO("TERMINATING\n");
}

#endif



