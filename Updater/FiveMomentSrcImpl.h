// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for five-moment source terms
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GK_FIVE_MOMENT_SRC_H
#define GK_FIVE_MOMENT_SRC_H

#include <stdint.h>

extern "C" {
    typedef struct {
        double charge, mass; /* Charge and mass */
    } FluidData_t;

    typedef struct {
        int8_t nFluids; /* Number of fluids */
        double epsilon0; /* Permittivity of free space */
        double chi_e, chi_m; /* Propagation speed factor for electric field error potential */
        int8_t gravityDir; /* Direction of gravity force */
        double gravity; /* Gravitational acceleration */
        bool hasStatic, hasPressure; /* Flag to indicate if there is: static EB field, pressure */
    } FiveMomentSrcData_t;

    /* Method to update fluids and flield using explicit RK3 method */
    void gkylFiveMomentSrcRk3(FiveMomentSrcData_t *sd, FluidData_t *fd, double dt, double **f, double *em);
}

#endif // GK_FIVE_MOMENT_SRC_H