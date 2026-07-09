#pragma once 

// #define EXTRA_DEBUG

#define SPIN 1.0    //<! value of spin (1/2 -> spin operators, 1 -> pauli matrices)
#define CONFIG 7    //<! on-site configuration (local hilbert space dimension)


// #define NO_POINT_SYM


#if defined(MY_MAC) // use only on personal device

    //<! Macro to control if main routines are set for symmetric or non-symmetric model
    #ifndef USE_SYMMETRIES
        // #define USE_SYMMETRIES
    #endif
    
    #undef USE_REAL_SECTORS
    //<! Macro to set element type to double for real momentum sectors
    #ifndef USE_REAL_SECTORS
        #define USE_REAL_SECTORS
    #endif

    //<! Macro to remove all other point symmetries
    #ifndef USE_ONLY_TRANSLATION
        // #define USE_ONLY_TRANSLATION
    #endif

#endif


#include "compiler.hpp"