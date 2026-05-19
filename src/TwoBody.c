#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "RK4.h"

#define SATELLITE_MASS 1000.0
#define ALTITUDE 300000.0
#define R0 (EARTH_RADIUS + ALTITUDE)

// Burn parameters
#define DELTA_V1 2000.0
#define DELTA_V2 1000.0
#define DELTA_V3 1500.0

#define DT 10.0
#define NUM_STEPS 12000    // enough to see all three burns

// Apply burn
// Prograde if dv > 0
// Retrograde if dv < 0
void apply_burn(Body *body, double dv) {
    double v = sqrt(body->vx*body->vx + body->vy*body->vy + body->vz*body->vz);
    double tx = body->vx / v;
    double ty = body->vy / v;
    double tz = body->vz / v;
    body->vx += tx * dv;
    body->vy += ty * dv;
    body->vz += tz * dv;
    if (dv > 0)
        printf("Prograde burn: dv = %.1f m/s\n", dv);
    else
        printf("Retrograde burn: dv = %.1f m/s\n", dv);
}

// Compute orbital parameters relative to Earth
void compute_relative_orbit(Body *earth, Body *sc, double mu, double *r_peri, double *r_apo) {
    double rx = sc->x - earth->x;
    double ry = sc->y - earth->y;
    double rz = sc->z - earth->z;
    double vx = sc->vx - earth->vx;
    double vy = sc->vy - earth->vy;
    double vz = sc->vz - earth->vz;
    double r = sqrt(rx*rx + ry*ry + rz*rz);
    double v = sqrt(vx*vx + vy*vy + vz*vz);
    double hx = ry*vz - rz*vy;
    double hy = rz*vx - rx*vz;
    double hz = rx*vy - ry*vx;
    double h = sqrt(hx*hx + hy*hy + hz*hz);
    double ex = (v*v*rx - (rx*vx + ry*vy + rz*vz)*vx) / mu - rx/r;
    double ey = (v*v*ry - (rx*vx + ry*vy + rz*vz)*vy) / mu - ry/r;
    double ez = (v*v*rz - (rx*vx + ry*vy + rz*vz)*vz) / mu - rz/r;
    double e = sqrt(ex*ex + ey*ey + ez*ez);
    double a = h*h / (mu * (1.0 - e*e));
    *r_peri = a * (1.0 - e);
    *r_apo  = a * (1.0 + e);
}



int main() {
    const double mu_earth = G * EARTH_MASS;

    // Two bodies: Earth (0) and satellite (1)
    Body bodies[2];

    // Initial conditions for circular LEO, barycenter at origin
    double v_rel_circ = sqrt(mu_earth / R0);
    double r_sat_rel = R0;
    double v_sat_rel_x = 0.0;
    double v_sat_rel_y = v_rel_circ;
    double v_sat_rel_z = 0.0;

    double M_total = EARTH_MASS + SATELLITE_MASS;
    double r_earth_x = -(SATELLITE_MASS / M_total) * r_sat_rel;
    double r_earth_y = 0.0;
    double r_earth_z = 0.0;
    double r_sat_x = (EARTH_MASS / M_total) * r_sat_rel;
    double r_sat_y = 0.0;
    double r_sat_z = 0.0;

    double v_earth_x = -(SATELLITE_MASS / M_total) * v_sat_rel_x;
    double v_earth_y = -(SATELLITE_MASS / M_total) * v_sat_rel_y;
    double v_earth_z = -(SATELLITE_MASS / M_total) * v_sat_rel_z;
    double v_sat_x = (EARTH_MASS / M_total) * v_sat_rel_x;
    double v_sat_y = (EARTH_MASS / M_total) * v_sat_rel_y;
    double v_sat_z = (EARTH_MASS / M_total) * v_sat_rel_z;

    bodies[0] = (Body){r_earth_x, r_earth_y, r_earth_z,
                       v_earth_x, v_earth_y, v_earth_z,
                       EARTH_MASS};
    bodies[1] = (Body){r_sat_x, r_sat_y, r_sat_z,
                       v_sat_x, v_sat_y, v_sat_z,
                       SATELLITE_MASS};

    // Burn 1: at initial perigee (prograde)
    printf("\nBurn 1 (perigee, prograde)\n");
    apply_burn(&bodies[1], DELTA_V1);
    double r_peri, r_apo;
    compute_relative_orbit(&bodies[0], &bodies[1], mu_earth, &r_peri, &r_apo);
    printf("After burn 1: perigee = %.1f km, apogee = %.1f km\n",
           r_peri/1000.0, r_apo/1000.0);

    // Open CSV
    FILE *fp = fopen("orbit.csv", "w");
    if (!fp) {
        fprintf(stderr, "Error: cannot create orbit.csv\n");
        return 1;
    }
    fprintf(fp, "time_s,earth_x_m,earth_y_m,earth_z_m,sat_x_m,sat_y_m,sat_z_m,sat_vx_m_s,sat_vy_m_s,sat_vz_m_s\n");

    // Simulation loop with apogee/perigee detection for burns 2 and 3
    double time = 0.0;
    int burn1_done = 1;       // already done
    int burn2_done = 0;
    int burn3_done = 0;
    double prev_r_rel = 0.0;
    int rising = 1;           // 1 = distance increasing, 0 = decreasing
    int perigee_count = 0;    // count perigees after burn 2

    for (int step = 0; step <= NUM_STEPS; step++) {
        // Write current state
        fprintf(fp, "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
                time,
                bodies[0].x, bodies[0].y, bodies[0].z,
                bodies[1].x, bodies[1].y, bodies[1].z,
                bodies[1].vx, bodies[1].vy, bodies[1].vz);

        if (step == NUM_STEPS) break;

        // Compute relative distance (Earth‑centered)
        double dx = bodies[1].x - bodies[0].x;
        double dy = bodies[1].y - bodies[0].y;
        double dz = bodies[1].z - bodies[0].z;
        double r_rel = sqrt(dx*dx + dy*dy + dz*dz);

        // Detect extrema
        if (step > 0) {
            if (rising && r_rel < prev_r_rel) {
                // Passed a maximum (apogee)
                if (burn1_done && !burn2_done) {
                    printf("\nBurn 2 (apogee, prograde) at t = %.1f s\n", time - DT);
                    apply_burn(&bodies[1], DELTA_V2);
                    compute_relative_orbit(&bodies[0], &bodies[1], mu_earth, &r_peri, &r_apo);
                    printf("After burn 2: perigee = %.1f km, apogee = %.1f km\n",
                           r_peri/1000.0, r_apo/1000.0);
                    burn2_done = 1;
                }
            } else if (!rising && r_rel > prev_r_rel) {
                // Passed a minimum (perigee)
                if (burn2_done && !burn3_done) {
                    perigee_count++;
                    if (perigee_count == 1) { 
                        printf("\nBurn 3 (perigee, retrograde) at t = %.1f s\n", time - DT);
                        apply_burn(&bodies[1], -DELTA_V3);
                        compute_relative_orbit(&bodies[0], &bodies[1], mu_earth, &r_peri, &r_apo);
                        printf("After burn 3: perigee = %.1f km, apogee = %.1f km\n",
                               r_peri/1000.0, r_apo/1000.0);
                        burn3_done = 1;
                    }
                }
            }
            rising = (r_rel > prev_r_rel);
        }
        prev_r_rel = r_rel;

        // RK4 step
        rk4_step(bodies, 2, DT);
        time += DT;

        // Check Earth impact
        if (r_rel <= EARTH_RADIUS) {
            printf("\nCrashed back to Earth at t = %.1f s\n", time);
            break;
        }

        if ((step+1) % 1000 == 0)
            printf("Step %d, time = %.1f s\n", step+1, time);
    }

    fclose(fp);
    printf("\nResults written to 'data/orbit.csv'\n");
    return 0;
}