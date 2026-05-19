#ifndef BURN_H
#define BURN_H

#include <math.h>
#include <stdio.h>

#include "RK4.h"


void apply_burn_vector(Body *body, double dv, double ux, double uy, double uz) {
    body->vx += ux * dv;
    body->vy += uy * dv;
    body->vz += uz * dv;
    printf("Burn: dv = %.1f m/s (%s)\n", dv, (dv > 0 ? "prograde" : "retrograde"));
}

// ------------------------------------------------------------------
// Apply tangential (prograde/retrograde) burn based on current velocity
// ------------------------------------------------------------------
void apply_tangential_burn(Body *body, double dv) {
    double v = sqrt(body->vx*body->vx + body->vy*body->vy + body->vz*body->vz);
    double ux = body->vx / v;
    double uy = body->vy / v;
    double uz = body->vz / v;
    apply_burn_vector(body, dv, ux, uy, uz);
}

#endif