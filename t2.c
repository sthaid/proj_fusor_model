// kinetic diameter;
// https://en.wikipedia.org/wiki/Kinetic_diameter
// 
//   d^2 = 1 / (pi * l * n)
//
//     d = kinetic diameter
//     l = mean free path
//     n = number density of particles

// XXX review
// XXX make physics.h
// XXX checkin
// XXX sdl graphics 
// XXX record the simulation, and playback mode
// XXX auto stop
// XXX improve print_state
// XXX combine state into one struct


// XXX time how long it takes
// XXX optimize by minimizing floating point operations
// XXX multithread
// XXX try float instead of double


#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>

#include <assert.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <sys/queue.h>

#include "util_misc.h"

//
// defines
//

#define NUM_SECTION               (1000)  // XXX name
#define NUM_SIM_PARTICLE          (1000000)  // XXX name  max vs num
#define CHAMBER_LENGTH            (10.0)         // m
#define CHAMBER_CROSS_SECTION     (.001 * .001)  // m^2
#define CHAMBER_VOLUME            (CHAMBER_LENGTH * CHAMBER_CROSS_SECTION)
#define SECTION_VOLUME            (CHAMBER_VOLUME / NUM_SECTION)
#define SECTION_LENGTH            (CHAMBER_LENGTH / NUM_SECTION)
#define NUM_REAL_PARTICLE         (5e20 * CHAMBER_VOLUME)   // XXX 5e20?

#define SECT_NUM(sp)              ((sp)->position / SECTION_LENGTH)

#define CHAMBER_LEFT_END_TEMP     200
#define CHAMBER_RIGHT_END_TEMP    400  
#define GUESS_TEMP                300

#define SECTION_POSITION(num)        ((num) * SECTION_LENGTH)

#define D2_AMU 4.03
#define AMU_TO_KG(x)    ((x) * 1.66054e-27)
#define k  1.38066e-23    // Boltzmann constant J/K
#define TEMPERATURE_TO_VELOCITY(T) (sqrt(3 * k * (T) / AMU_TO_KG(D2_AMU)))
#define VELOCITY_TO_TEMPERATURE(v) (AMU_TO_KG(D2_AMU) * (v) * (v) / (3. * k))

//
// typedefs
//

typedef struct {
    struct {
        int64_t num_sim_particle;  // number of simulated particles in this section
        double  temperature;       // average particle temperature in this section
    } section[NUM_SECTION];
} state_t;

typedef struct {
    double position;
    double velocity;
} sim_particle_t;

// 
// variables
//

state_t        state;
int32_t        num_sim_particle;
sim_particle_t sim_particle[NUM_SIM_PARTICLE];

//
// prototypes
//

void test(void);
void init(void);
void simulate(void);
void print_state(double t);

// -----------------  MAIN  ----------------------------------------------------------------------------

int main(int argc, char **argv)
{
#if 1
    double number_density, d, mean_free_path;
    number_density = 4.92234e+20;  // NUM_D2_MOLECULES/CU_M
    d = 289e-12;   // meters
                   // 289 pm (picometer)  (for H2 molecule)
                   // a picometer is 1/1000 of a nanometer
    mean_free_path = 1 / (M_PI * d * d * number_density);
    printf("mean_free_path = %lf meters\n", mean_free_path);

    printf("VELOCITY AT %d Kelvin %lg m/s\n",
        GUESS_TEMP, TEMPERATURE_TO_VELOCITY(GUESS_TEMP));

    printf("max DeltaT = %lg  secs\n",
        mean_free_path / TEMPERATURE_TO_VELOCITY(GUESS_TEMP) / 10);
    
    printf("\n");
#endif

    INFO("HELLO\n");
    test();

    return 0;
}

// -----------------  SIMULATION TEST  -----------------------------------------------------------------

void test(void)
{
    uint64_t start_us, duration_us;

    // init 
    init();

    // simulate the particles
    start_us = microsec_timer();
    simulate();
    duration_us = microsec_timer() - start_us;
    INFO("DURATION %0.3lf secs\n", (double)duration_us/1000000);
}

void init(void)
{
    int32_t i, num_sim_particle_in_section;

    num_sim_particle_in_section = NUM_SIM_PARTICLE / NUM_SECTION;
    for (i = 0; i < NUM_SECTION; i++) {
        state.section[i].num_sim_particle = num_sim_particle_in_section;
        state.section[i].temperature      = GUESS_TEMP;
    }

    num_sim_particle = 0;
    for (i = 0; i < NUM_SECTION; i++) {
        double velocity, spacing, position;
        int32_t j;

        velocity = TEMPERATURE_TO_VELOCITY(state.section[i].temperature);
        spacing = SECTION_LENGTH / num_sim_particle_in_section;
        position = SECTION_POSITION(i) + spacing / 2;
        for (j = 0; j < num_sim_particle_in_section; j++) {   
            assert(num_sim_particle < NUM_SIM_PARTICLE);
            sim_particle[num_sim_particle].position = position;
            sim_particle[num_sim_particle].velocity = ((j & 1) ? velocity : -velocity);
            position += spacing;
            num_sim_particle++;
        }
    }

    INFO("num_sim_particle = %d\n\n", num_sim_particle);
    assert(num_sim_particle <= NUM_SIM_PARTICLE);
    assert(num_sim_particle > NUM_SIM_PARTICLE*0.99);   
}

void simulate(void)
{
    int32_t i;
    //int32_t coll_cnt;
    double t;

    static struct {
        int64_t off_left;
        int64_t off_right;
        int64_t switch_sects;
        int64_t check_for_col;
        int64_t got_col;
    } stats;

    #define DELTA_T  1e-6
    #define MAX_T    0.000015

    for (t = 0; t < MAX_T; t += DELTA_T) {
        //coll_cnt = 0;
        for (i = 0; i < num_sim_particle; i++) {
            sim_particle_t * sp = &sim_particle[i];
            sim_particle_t   sp_orig = *sp;

            int32_t orig_sect_num = SECT_NUM(&sp_orig);
            assert(orig_sect_num >= 0 && orig_sect_num < NUM_SECTION);

            // update position
            sp->position = sp->position + sp->velocity * DELTA_T;

            // if the particle is outside the chamber on the left then
            //   set poistion to 0
            //   set temperature to that of the left chamber wall
            //   set positive velocity based on temperature
            // endif
            if (sp->position < 0) {
                stats.off_left++;
                double temp_sum, sp_delta_temp;
                sp_delta_temp = CHAMBER_LEFT_END_TEMP - VELOCITY_TO_TEMPERATURE(sp->velocity);
                sp->position = 0;
                sp->velocity = TEMPERATURE_TO_VELOCITY(CHAMBER_LEFT_END_TEMP);
                temp_sum = state.section[orig_sect_num].num_sim_particle *
                           state.section[orig_sect_num].temperature;
                state.section[orig_sect_num].temperature =
                    (temp_sum + sp_delta_temp) /
                    state.section[orig_sect_num].num_sim_particle;
            }

            // if the particu is outside the chamber on the right then
            //   set poistion to CHAMBER_LENGTH
            //   set temperature to that of the right chamber wall
            //   set negative velocity based on temperature
            // endif
            if (sp->position >= CHAMBER_LENGTH) {
                stats.off_right++;
                double temp_sum, sp_delta_temp;
                sp_delta_temp = CHAMBER_RIGHT_END_TEMP - VELOCITY_TO_TEMPERATURE(sp->velocity);
                sp->position = CHAMBER_LENGTH - 1e-9;
                sp->velocity = -TEMPERATURE_TO_VELOCITY(CHAMBER_RIGHT_END_TEMP);
                temp_sum = state.section[orig_sect_num].num_sim_particle *
                           state.section[orig_sect_num].temperature;
                state.section[orig_sect_num].temperature =
                    (temp_sum + sp_delta_temp) /
                    state.section[orig_sect_num].num_sim_particle;
            }


            int32_t new_sect_num = SECT_NUM(sp);
            assert(new_sect_num >= 0 && new_sect_num < NUM_SECTION);

            // if the particle has moved to a new section then update
            // the state of the old and new section
            if (new_sect_num != orig_sect_num) {
                stats.switch_sects++;
                double temp_sum;
                double sp_temperature = VELOCITY_TO_TEMPERATURE(sp->velocity);

// XXX frequent
// - need sp->temperature
                temp_sum = state.section[orig_sect_num].num_sim_particle *
                           state.section[orig_sect_num].temperature;
                state.section[orig_sect_num].temperature =
                    (temp_sum - sp_temperature) /
                    (state.section[orig_sect_num].num_sim_particle - 1);

                temp_sum = state.section[new_sect_num].num_sim_particle *
                           state.section[new_sect_num].temperature;
                state.section[new_sect_num].temperature =
                    (temp_sum + sp_temperature) /
                    (state.section[new_sect_num].num_sim_particle + 1);

                state.section[orig_sect_num].num_sim_particle--;
                state.section[new_sect_num].num_sim_particle++;
            }


            // based on the probability of a collision with another particle determine
            // whether a collision occurred; if so then update velocity using conservation
            // of momentum and the velocity of the particle and the velocity of particles
            // in this section
            // https://en.wikipedia.org/wiki/Kinetic_diameter
// XXX frequent
// - need optimized exp

            stats.check_for_col++;
            double d, mean_free_path, average_num_events_per_sec, average_num_events_per_interval;
            double number_density, probability_of_0_collisions_in_interval;
            bool collision_occurred;

            number_density = state.section[new_sect_num].num_sim_particle *
                             (NUM_REAL_PARTICLE / num_sim_particle) /
                             SECTION_VOLUME;
            d = 289e-12;   // kinetic diameter of hydrogen in meters for H2
            mean_free_path = 1. / (M_PI * d * d * number_density);
            average_num_events_per_sec = fabs(sp->velocity) / mean_free_path;
            average_num_events_per_interval = average_num_events_per_sec * DELTA_T;
            probability_of_0_collisions_in_interval =   exp(-average_num_events_per_interval);
            collision_occurred = ((double)random()/RAND_MAX > probability_of_0_collisions_in_interval);

#if 0
            printf("sp[%d]: position=%lf velocity=%lf temp=%lf\n",
                   i, sp->position, sp->velocity, VELOCITY_TO_TEMPERATURE(sp->velocity));
            printf("  new_sect_num=%d  section_num_density=%lg\n",
                   new_sect_num, number_density);
            printf("  mfp=%lf  avg_col/sec=%lf  avg_col_per_intvl=%lf\n",
                   mean_free_path, average_num_events_per_sec, average_num_events_per_interval);
            printf("  probability_of_0_col_in_intvl = %lf\n",
                   probability_of_0_collisions_in_interval);
#endif
            if (collision_occurred) {
                stats.got_col++;
#if 0
                coll_cnt++;
                printf("*** COLLISION i=%d coll_cnt=%d fract=%lf ***\n",
                       i, coll_cnt, (double)coll_cnt / (i+1));
#endif
                // assume particles are moving towards each other, and
                // they exchange velocities due to conservation of momentum
                double vs;
// XXX frequent
// need section velocity
                vs = TEMPERATURE_TO_VELOCITY(state.section[new_sect_num].temperature);
                if (sp->velocity > 0) {
                    sp->velocity = -vs;
                } else {
                    sp->velocity = vs;
                }
            }
        }

#if 0
        static int64_t print_stats_count;
        if ((print_stats_count % 10) == 0) {
            print_state(t);
        }
        print_stats_count++;
#endif
    }

    printf("%8ld %8ld %8ld %8ld %8ld\n",
        stats.off_left,
        stats.off_right,
        stats.switch_sects,
        stats.check_for_col,
        stats.got_col);
    printf("\n");

    print_state(t);
}

void print_state(double time)
{
    int32_t i;

#if 0
    printf("%s\n", s);
    for (i = 0; i < NUM_SECTION; i++) {
        printf("%s: sect %d:  num_sim_particle=%ld  temperature=%lf\n", 
           __func__, i, state.section[i].num_sim_particle, state.section[i].temperature);
    }
    printf("\n");
#endif
#if 0
    printf("%20s: %10.3f %10.3f %10.3f %10.3f %10.3f    NUM_SECTION=%d\n", 
        s,
        state.section[0].temperature,  
        state.section[NUM_SECTION*1/4].temperature,  
        state.section[NUM_SECTION*2/4].temperature,  
        state.section[NUM_SECTION*3/4].temperature,  
        state.section[NUM_SECTION-1].temperature, 
        NUM_SECTION);
#endif

    printf("%8.6lf ", time);
    for (i = 0; i < NUM_SECTION; i++) {
        double t = state.section[i].temperature;
        printf("%8.3f ", t);
        if (fabs(t-300) < .001) break;
    }
    printf("   I=%d\n", i);
    printf("\n");

}
