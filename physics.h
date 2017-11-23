#ifndef __PHYSICS_H__
#define __PHYSICS_H__

// conversion macros
#define AMU_TO_KG(amu)          ((amu) * 1.66054e-27)            // mass 
#define MTORR_TO_PASCAL(mtorr)  ((mtorr) * 0.13332237)           // pressure
#define UTORR_TO_PASCAL(utorr)  ((utorr) * 0.00013332237)
#define F_TO_C(tf)              (((tf) - 32.0) * (5.0 / 9.0))    // temperature
#define C_TO_K(tc)              ((tc) + 273.15)
#define F_TO_K(tf)              (C_TO_K(F_TO_C(tf)))
#define M_TO_IN(m)              ((m) * 39.3701)                  // length
#define M_TO_NM(m)              ((m)  * 1000000000)
#define NM_TO_M(nm)             ((nm) / 1000000000)
#define M_TO_MM(m)              ((m)  * 1000)
#define MM_TO_M(mm)             ((mm) / 1000)
#define MM_TO_NM(mm)            ((mm) * 1000000)
#define NM_TO_MM(nm)            ((nm) / 1000000)
#define MPERS_TO_NMPERNS(v)     (v)                             // velocity
#define JOULES_TO_EV(j)         ((j) * 6.242e+18)               // energy
#define EV_TO_JOULES(ev)        ((ev) * 1.60218e-19)

// temperature 
#define ROOM_TEMPERATURE_K 293.0

// deuterium 
#define D_AMU   2.000  // xxx check this
#define D_MASS  AMU_TO_KG(D_AMU)

// hydrogen
#define H_IONIZATION_ENERGY_EV  13.5984
#define H2_KINETIC_DIAMETER     289e-12 // meters

// electron 
#define ELECTRON_CHARGE -1.60217662e-19   // C
#define ELECTRON_MASS    9.10938356e-31   // kg

// proton
#define PROTON_CHARGE   -ELECTRON_CHARGE  // C
#define PROTON_MASS     1.6726219e-27     // kg

// kinetic temperature
// reference: http://hyperphysics.phy-astr.gsu.edu/hbase/Kinetic/kintem.html
#define K 1.38066e-23  // Boltzmann constant Joules/Kelvin
#define TEMPERATURE_TO_VELOCITY(t,m) (sqrt((t) * (3. * K / (m))))
#define VELOCITY_TO_TEMPERATURE(v,m) (((m) / (3. * K)) * (v) * (v))

// kinetic energy
#define KINETIC_ENERGY(m,v)              (0.5 * (m) * (v) * (v))
#define KINETIC_ENERGY_TO_VELOCITY(ke,m) (sqrt(2. * (ke) / (m)))

// ideal gas law
// http://hyperphysics.phy-astr.gsu.edu/hbase/Kinetic/idegas.html
#define NUMBER_OF_PARTICLES(p,v,t) ((p) * (v) / (K * (t)))
#define NUMBER_DENSITY(p,t)        ((p) / (K * (t)))

// Coulomb's Law
// http://hyperphysics.phy-astr.gsu.edu/hbase/electric/elefor.html
//    F = KE q1 * q2 / r^2
#define KE  8.987552e9  // Coulomb's Constant  N m^2 C^-2
#define COULOMB_FORCE(q1,q2,r)  ((KE * (q1) * (q2)) / ((r) * (r)))  // repulsive

// Electric Potential
// http://physics.bu.edu/~duffy/PY106/Potential.html
//  V = KE Q / r
//  Q = V * r / KE
#define CHARGE_TO_POTENTIAL(q,r) (KE * (q) / (r))
#define POTENTIAL_TO_CHARGE(v,r)  ((v) * (r) / KE)

// Mean Free Path, for stationary target
// http://hyperphysics.phy-astr.gsu.edu/hbase/Kinetic/menfre.html
// https://en.wikipedia.org/wiki/Mean_free_path
// https://en.wikipedia.org/wiki/Kinetic_diameter
//                             1
// mean_free_path = -----------------------------
//                  cross_section * number_density
#define CROSS_SECTION(d)       (M_PI * (d) * (d))     // d = kinetic diameter
#define MEAN_FREE_PATH(cs,nd)  (1.0 / ((cs) * (nd)))

#endif
