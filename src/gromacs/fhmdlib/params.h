#ifndef FHMD_PARAMS_H_
#define FHMD_PARAMS_H_

#define FHMD_VERSION            1.00            /* FHMD model version */

#define FHMD_MAX_LENGTH         1000.0          /* Maximum length scale [nm] -- for control purpose only */

#define FHMD_kB                 0.00831451      /* Boltzmann constant [kJ/(mol*K)] */

//#define FHMD_DEBUG                            /* Write debug information */
//#define FHMD_TECPLOT                          /* Write data to Tecplot (create 'tecplot' dir manually) */
//#define FHMD_DEBUG_GRID                       /* Print FH grid coordinates to the screen */
//#define FHMD_DEBUG_FH                         /* Print FH debug information */
//#define FHMD_DEBUG_COM                        /* Print protein COM coordinates */
//#define FHMD_DEBUG_INTERPOL                   /* Print interpolated values for debugging */

enum FHMD_SCHEME {Pure_MD, One_Way, Two_Way};   /* List of schemes */

enum FHMD_EOS    {eos_argon, eos_spce};         /* Equation of state enumeration */

enum FHMD_S      {constant_S, fixed_sphere, moving_sphere};

enum FHMD_CELL   {FH_zone, hybrid_zone, boundary};

/* EOS coefficients for RIGID SPC/E water */
#define FHMD_EOS_SPCE_MU        409.496
#define FHMD_EOS_SPCE_KAPPA     933.41
#define FHMD_EOS_SPCE_A         0.010102137
#define FHMD_EOS_SPCE_B        -10.133069
#define FHMD_EOS_SPCE_C         2428.9203

/* EOS coefficients for Liquid Argon (OPLS-AA opls_097 force field) */
//#define FHMD_MASS_AR            39.948          /* Mass of Argon atom, amu */
//#define FHMD_SIGMA_AR           3.40100e-01     /* L-J Sigma parameter, nm */
//#define FHMD_EPSILON_AR         9.78638e-01     /* L-J Epsilon parameter, kJ/mol */
//#define FHMD_EOS_ARGON_MU       (1.011213*(FHMD_MASS_AR/(FHMD_SIGMA_AR*FHMD_TIME_AR)))
//#define FHMD_EOS_ARGON_KAPPA    (0.336763*(FHMD_MASS_AR/(FHMD_SIGMA_AR*FHMD_TIME_AR)))
//#define FHMD_EOS_ARGON_A        1.123600
//#define FHMD_EOS_ARGON_B        435.5957
//#define FHMD_EOS_ARGON_C        4.14560e+07
//#define FHMD_EOS_ARGON_D        (2.39206e-08*(FHMD_MASS_AR/(FHMD_SIGMA_AR*FHMD_TIME_AR*FHMD_TIME_AR)))
//#define FHMD_EOS_ARGON_E        (1641.000*((FHMD_SIGMA_AR*FHMD_SIGMA_AR*FHMD_SIGMA_AR)/FHMD_MASS_AR))

// from Voulgarakis:
#define FHMD_MASS_AR            39.948          /* Mass of Argon atom, amu */
#define FHMD_SIGMA_AR           3.40500e-01     /* L-J Sigma parameter, nm */
#define FHMD_EPSILON_AR         9.95792e-01     /* L-J Epsilon parameter, kJ/mol */
#define FHMD_TIME_AR            (FHMD_SIGMA_AR*sqrt(FHMD_MASS_AR/FHMD_EPSILON_AR))
#define FHMD_EOS_ARGON_MU       (1.011213*(FHMD_MASS_AR/(FHMD_SIGMA_AR*FHMD_TIME_AR)))
#define FHMD_EOS_ARGON_KAPPA    (0.336763*(FHMD_MASS_AR/(FHMD_SIGMA_AR*FHMD_TIME_AR)))
/* p in gmx units, exp function */
#define FHMD_EOS_ARGON_A       -34.85802
#define FHMD_EOS_ARGON_B        0.01932
#define FHMD_EOS_ARGON_C        110.76002
#define FHMD_EOS_ARGON_D        0.0
#define FHMD_EOS_ARGON_E        0.0
/* from close to Voulgarakis (Berendsen thermostat and 0.9-1 nm LJ switching):
#define FHMD_EOS_ARGON_A       -34.57909
#define FHMD_EOS_ARGON_B        0.02944
#define FHMD_EOS_ARGON_C        116.06787
#define FHMD_EOS_ARGON_D        0.0
#define FHMD_EOS_ARGON_E        0.0 
*/
/* Color codes */
#ifndef _MSC_VER
    #define RESET_COLOR         "\e[m"
    #define MAKE_RED            "\e[91m"
    #define MAKE_GREEN          "\e[92m"
    #define MAKE_YELLOW         "\e[93m"
    #define MAKE_BLUE           "\e[94m"
    #define MAKE_PURPLE         "\e[95m"
    #define MAKE_LIGHT_BLUE     "\e[96m"
#else
    #define RESET_COLOR
    #define MAKE_RED
    #define MAKE_GREEN
    #define MAKE_YELLOW
    #define MAKE_BLUE
    #define MAKE_PURPLE
    #define MAKE_LIGHT_BLUE
#endif

#endif /* FHMD_PARAMS_H_ */

