// XXX review all this

//
// CONVERSIONS
//

// length
#define IN_TO_M(len)            ((len) * 0.0254)
#define M_TO_IN(len)            ((len) / 0.0254)

// volume
#define IN3_TO_M3(vol)          ((vol) * 1.63871e-5)
#define M3_TO_C3(vol)           ((vol) * 1e6)

// temperature
#define F_TO_C(t)               (((t) - 32.) * (5. / 9.))
#define F_TO_K(t)               (F_TO_C(t) + 273.15)

// mass
#define AMU_TO_KG(amu)          ((amu) * 1.66054e-27)

// pressure
#define MTORR_TO_PASCAL(mtorr)  ((mtorr) * 0.133322) 

//
// CONSTANTS
//

#define k      1.38066e-23         // Boltzmann constant Joules/Kelvin
#define D2_AMU 4.03                // Deuterium molecule mass
#define D2_KG  AMU_TO_KG(D2_AMU)

//
// KINETIC TEMPERATURE
//
// http://hyperphysics.phy-astr.gsu.edu/hbase/Kinetic/kintem.html
//

#define TEMPERATURE_TO_VELOCITY(t,m) (sqrt((t) * (3. * k / (m))))
#define VELOCITY_TO_TEMPERATURE(v,m) (((m) / (3. * k)) * (v) * (v))

//
// IDEAL GAS LAW
//
// http://hyperphysics.phy-astr.gsu.edu/hbase/Kinetic/idegas.html
//

#define NUMBER_OF_MOLECULES(p,v,t)        ((p) * (v) / (k * (t)))
#define NUMBER_DENSITY_OF_MOLECULES(p,t)  ((p) / (k * (t)))

//
// SI UNITS
// 

// XXX tbd

//
// OTHER REFERENCES
//

// Deuterium
// http://www.sigmaaldrich.com/catalog/product/aldrich/368407?lang=en&region=US


// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#if 0
// Poisson Distribution
// https://en.wikipedia.org/wiki/Poisson_distribution

// Mean Free Path
// http://hyperphysics.phy-astr.gsu.edu/hbase/Kinetic/menfre.html

// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// kinetic diameter;
// https://en.wikipedia.org/wiki/Kinetic_diameter
// 
//   d^2 = 1 / (pi * l * n)
//
//     d = kinetic diameter
//     l = mean free path
//     n = number density of particles






#define R    8.3145   // universal gas constant  J/mol K
#define NA   6.0221E23   // avogadro's number

#define k  1.38066e-23    // Boltzmann constant J/K

#define AMU_TO_KG(x)    ((x) * 1.66054e-27)             // xxx

// Isotopes of Hydrogen
//https://www.boundless.com/chemistry/textbooks/boundless-chemistry-textbook/nonmetallic-elements-21/hydrogen-148/isotopes-of-hydrogen-573-3644/
#define D2_AMU 4.03

#define ELECTRON_CHARGE  1.60217662e-19   // Coulombs or C
#define ELECTRON_MASS    9.10938356e-31   // kilograms  or kg

#define GRID_VOLTAGE 30000.0
#define GRID_RADIUS  IN_TO_M(0.75)   
#define CHAMBER_RADIUS IN_TO_M(3)

//  1 Amp = 1 Coulomb / Sec


// http://physics.bu.edu/~duffy/PY106/Potential.html
//  V = ke Q / r

// Coulomb's Law
//    F = ke q1 * q2 / r^2
//      (ke = 8.99×10^9 N m2 C−2)
//     8.987552e9  N m^2 C^-2

#define ke  8.987552e9  // Coulomb's Constant  N m^2 C^-2

#define NANOSECOND  1e-9

#define JOULES_TO_EV(x)    ((x) * 6.242e+18)

#define h   6.62607004e-34 // m2 kg / s   Planck constant
#define c   299792458.0    //  m / s      Speed of light

#define M_TO_NM(x)   ((x) * 1e9)

#define ELECTRONS_PER_COULOMB  (1. / ELECTRON_CHARGE)
#define COULOMB_PER_ELECTRON   (ELECTRON_CHARGE)

#endif
