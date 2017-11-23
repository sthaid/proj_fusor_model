#if 0
XXX redo the results, and comment

[HOME haid@home proj_fusor_model]$ ./mfp_h2  15 300
pressure       = 15.000000 mtorr
pressure       = 1.999836 pascal
temperature    = 300.000000 kelvin  3.000000e+02
velocity       = 1362.659763 m/s  1.362660e+03
kinetic energy = 0.038781 ev
cs             = 2.623890e-19 m^2
nd             = 4.943561e+20 m^-3
mfp            = 0.007709 m
time_mfp       = 5.657533 us

[HOME haid@home proj_fusor_model]$ ./mfp_h2  15 400000000
pressure       = 15.000000 mtorr
pressure       = 1.999836 pascal
temperature    = 400000000.000000 kelvin  4.000000e+08
velocity       = 1573463.962400 m/s  1.573464e+06
kinetic energy = 51708.478320 ev
cs             = 2.623890e-19 m^2
nd             = 4.943561e+20 m^-3
mfp            = 0.007709 m
time_mfp       = 0.004900 us

[HOME haid@home proj_fusor_model]$ ./mfp_h2  1000  400000000
pressure       = 1000.000000 mtorr
pressure       = 133.322370 pascal
temperature    = 400000000.000000 kelvin  4.000000e+08
velocity       = 1573463.962400 m/s  1.573464e+06
kinetic energy = 51708.478320 ev
cs             = 2.623890e-19 m^2
nd             = 3.295708e+22 m^-3
mfp            = 0.000116 m
time_mfp       = 0.000073 us
#endif

#include <stdio.h>
#include <math.h>
#include "physics.h"

int main(int argc, char **argv) 
{
    double pressure_mtorr, temp_kelvin;
    double pressure_pascal, cs, nd, mfp, velocity, time_mfp;
    double ke_joules, ke_ev;

    if (argc != 3 || 
        sscanf(argv[1], "%lf", &pressure_mtorr) != 1 ||
        sscanf(argv[2], "%lf", &temp_kelvin) != 1)
    {
        printf("usage: mfp_h2 <mtorr> <temp kelvin>\n");
        return 1;
    }

    // XXX mixes D and H2
    pressure_pascal = MTORR_TO_PASCAL(pressure_mtorr);
    velocity = TEMPERATURE_TO_VELOCITY(temp_kelvin, D_MASS);
    ke_joules = KINETIC_ENERGY(D_MASS, velocity);
    ke_ev = JOULES_TO_EV(ke_joules);
    cs = CROSS_SECTION(H2_KINETIC_DIAMETER);
    nd = NUMBER_DENSITY(pressure_pascal,ROOM_TEMPERATURE_K);
    mfp = MEAN_FREE_PATH(cs,nd);
    time_mfp = mfp / velocity;
   
    printf("pressure       = %f mtorr\n", pressure_mtorr);
    printf("pressure       = %f pascal\n", pressure_pascal);
    printf("temperature    = %f kelvin  %e\n", temp_kelvin, temp_kelvin);
    printf("velocity       = %f m/s  %e\n", velocity, velocity);
    printf("kinetic energy = %f ev\n", ke_ev);
    printf("cs             = %e m^2\n", cs);
    printf("nd             = %e m^-3\n", nd);
    printf("mfp            = %f m\n", mfp);
    printf("time_mfp       = %f us\n", time_mfp * 1e6);         

    return 0;
}
