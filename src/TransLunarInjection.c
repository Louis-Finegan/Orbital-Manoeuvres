#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "RK4.h"
#include "Burn.h"


/* Satellite initial orbit (circular LEO) */
#define ALTITUDE     300000.0
#define R0           (EARTH_RADIUS + ALTITUDE)

/* Moon orbit parameters */
#define MOON_DIST    384400000.0   /* average Earth‑Moon distance (m) */
#define MOON_SPEED   1022.0        /* orbital speed of Moon around Earth (m/s) */

/* TLI burn – tangential prograde at perigee */
#define TLI_DV       3150.0        /* m/s – adjust for closer flyby */

/* Lunar orbit insertion (LOI) target altitude (m) */
#define DESIRED_LUNAR_ALTITUDE 200000.0   /* 100 km above Moon's surface */

/* Moon initial phase angle (degrees) – adjust for intercept */
#define MOON_PHASE_DEG 133.83

/* Simulation parameters */
#define DT           10.0
#define NUM_STEPS    864000



// Compute orbit relative to a central body
void compute_orbit_relative_to(Body *center, Body *sat, double mu_center,
                               double *r_peri, double *r_apo) {
    double rx = sat->x - center->x;
    double ry = sat->y - center->y;
    double rz = sat->z - center->z;
    double vx = sat->vx - center->vx;
    double vy = sat->vy - center->vy;
    double vz = sat->vz - center->vz;
    double r = sqrt(rx*rx + ry*ry + rz*rz);
    double v = sqrt(vx*vx + vy*vy + vz*vz);
    double hx = ry*vz - rz*vy;
    double hy = rz*vx - rx*vz;
    double hz = rx*vy - ry*vx;
    double h = sqrt(hx*hx + hy*hy + hz*hz);
    double ex = (v*v*rx - (rx*vx + ry*vy + rz*vz)*vx) / mu_center - rx/r;
    double ey = (v*v*ry - (rx*vx + ry*vy + rz*vz)*vy) / mu_center - ry/r;
    double ez = (v*v*rz - (rx*vx + ry*vy + rz*vz)*vz) / mu_center - rz/r;
    double e = sqrt(ex*ex + ey*ey + ez*ez);
    double a = h*h / (mu_center * (1.0 - e*e));
    *r_peri = a * (1.0 - e);
    *r_apo  = a * (1.0 + e);
}

// ------------------------------------------------------------------
int main() {
    const double mu_earth = G * EARTH_MASS;
    const double mu_moon  = G * MOON_MASS;
    const double moon_phase_rad = MOON_PHASE_DEG * M_PI / 180.0;
    const double target_orbit_radius = MOON_RADIUS + DESIRED_LUNAR_ALTITUDE;

    // Three bodies: Earth (0), Moon (1) and Satellite (2)
    Body bodies[3];

    // Barycentric initial conditions for Earth and Moon
    double M_total = EARTH_MASS + MOON_MASS;
    double r_earth_bary = -(MOON_MASS / M_total) * MOON_DIST;
    double r_moon_bary  =  (EARTH_MASS / M_total) * MOON_DIST;

    // Moon initial position
    double moon_x0 = r_moon_bary * cos(moon_phase_rad);
    double moon_y0 = r_moon_bary * sin(moon_phase_rad);
    double moon_z0 = 0.0;
    double moon_vx0 = -MOON_SPEED * sin(moon_phase_rad);
    double moon_vy0 =  MOON_SPEED * cos(moon_phase_rad);
    double moon_vz0 = 0.0;

    // Earth initial position
    double earth_x0 = r_earth_bary * cos(moon_phase_rad);
    double earth_y0 = r_earth_bary * sin(moon_phase_rad);
    double earth_z0 = 0.0;
    double earth_vx0 = -(MOON_MASS / EARTH_MASS) * moon_vx0;
    double earth_vy0 = -(MOON_MASS / EARTH_MASS) * moon_vy0;
    double earth_vz0 = 0.0;

    bodies[0] = (Body){earth_x0, earth_y0, earth_z0,
                       earth_vx0, earth_vy0, earth_vz0,
                       EARTH_MASS};
    bodies[1] = (Body){moon_x0, moon_y0, moon_z0,
                       moon_vx0, moon_vy0, moon_vz0,
                       MOON_MASS};

    // Satellite initial condition: circular LEO around Earth
    double sat_rel_x = R0;
    double sat_rel_y = 0.0;
    double sat_rel_z = 0.0;
    double sat_rel_vx = 0.0;
    double sat_rel_vy = sqrt(mu_earth / R0);
    double sat_rel_vz = 0.0;

    double sat_x0 = earth_x0 + sat_rel_x;
    double sat_y0 = earth_y0 + sat_rel_y;
    double sat_z0 = earth_z0 + sat_rel_z;
    double sat_vx0 = earth_vx0 + sat_rel_vx;
    double sat_vy0 = earth_vy0 + sat_rel_vy;
    double sat_vz0 = earth_vz0 + sat_rel_vz;

    bodies[2] = (Body){sat_x0, sat_y0, sat_z0,
                       sat_vx0, sat_vy0, sat_vz0,
                       SATELLITE_MASS};

    // TLI burn (prograde)
    printf("\nTrans-Lunar Injection (TLI) burn\n");
    apply_tangential_burn(&bodies[2], TLI_DV);

    double r_peri, r_apo;

    compute_orbit_relative_to(&bodies[0], &bodies[2], mu_earth, &r_peri, &r_apo);
    printf("Earth-relative orbit after TLI:\n");
    printf("Perigee = %.1f km, Apogee = %.1f km\n", r_peri/1000.0, r_apo/1000.0);

    if (r_apo > MOON_DIST)
        printf("Apogee exceeds Moon's orbital radius! Lunar intercept possible.\n");
    else
        printf("Warning: apogee below Moon's orbit! May not reach Moon.\n");

    // Open CSV file
    FILE *fp = fopen("data/lunar_transfer.csv", "w");
    if (!fp) {
        fprintf(stderr, "Error: cannot create lunar_transfer.csv\n");
        return 1;
    }
    fprintf(fp, "time_s,earth_x,earth_y,earth_z,moon_x,moon_y,moon_z,sat_x,sat_y,sat_z\n");

    // --- Simulation loop with LOI burn at perilune ---
    double time = 0.0;
    int loi_done = 0;
    double closest_approach = 1e30;
    double prev_r = 1e30;
    double r_at_burn = 0.0, v_rel_at_burn = 0.0;
    double burn_dv = 0.0;

    for (int step = 0; step <= NUM_STEPS; step++) {
        // Write current state
        fprintf(fp, "%.3f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f\n",
                time,
                bodies[0].x, bodies[0].y, bodies[0].z,
                bodies[1].x, bodies[1].y, bodies[1].z,
                bodies[2].x, bodies[2].y, bodies[2].z);

        if (step == NUM_STEPS) break;

        // Relative position and velocity to Moon
        double dx = bodies[2].x - bodies[1].x;
        double dy = bodies[2].y - bodies[1].y;
        double dz = bodies[2].z - bodies[1].z;
        double dvx = bodies[2].vx - bodies[1].vx;
        double dvy = bodies[2].vy - bodies[1].vy;
        double dvz = bodies[2].vz - bodies[1].vz;
        double r_sat_moon = sqrt(dx*dx + dy*dy + dz*dz);
        double v_rel_moon = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);

        // Update closest approach info before potential burn
        if (r_sat_moon < closest_approach) {
            closest_approach = r_sat_moon;
            r_at_burn = r_sat_moon;
            v_rel_at_burn = v_rel_moon;
        }

        // Detect lunar perilune
        if (!loi_done && step > 0 && r_sat_moon > prev_r && prev_r < 1e29) {
            // Apply LOI burn
            printf("Lunar closest approach detected at t = %.1f s (%.2f days)\n",
                   time - DT, (time - DT)/86400.0);
            printf("Closest approach distance = %.1f km\n", closest_approach/1000.0);
            printf("Relative speed at that point = %.1f m/s\n", v_rel_at_burn);

            // Compute required circular speed at the target orbit radius
            double v_circ_target = sqrt(mu_moon / target_orbit_radius);
            double v_circ_at_r = sqrt(mu_moon / r_at_burn);

            burn_dv = v_rel_at_burn - v_circ_at_r;
            if (burn_dv > 0) {
                printf("Required retrograde burn = %.1f m/s to circularise at %.1f km altitude\n", burn_dv, DESIRED_LUNAR_ALTITUDE/1000.0);

                // Burn direction (retrograde)
                double inv_vx = -dvx / v_rel_at_burn;
                double inv_vy = -dvy / v_rel_at_burn;
                double inv_vz = -dvz / v_rel_at_burn;
                apply_burn_vector(&bodies[2], burn_dv, inv_vx, inv_vy, inv_vz);
                loi_done = 1;
                // After burn, recompute relative speed for info
                dvx = bodies[2].vx - bodies[1].vx;
                dvy = bodies[2].vy - bodies[1].vy;
                dvz = bodies[2].vz - bodies[1].vz;
                v_rel_moon = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
                printf("New relative speed = %.1f m/s\n", v_rel_moon);
            } else {
                printf("No burn needed! Already at or below circular speed.\n");
                loi_done = 1;
            }
        }
        prev_r = r_sat_moon;

        // Moon impact check
        if (r_sat_moon <= MOON_RADIUS) {
            printf("Lunar impact at t = %.1f s (%.2f days)\n", time, time/86400.0);
            break;
        }

        // Earth impact check
        double dxe = bodies[2].x - bodies[0].x;
        double dye = bodies[2].y - bodies[0].y;
        double dze = bodies[2].z - bodies[0].z;
        double r_sat_earth = sqrt(dxe*dxe + dye*dye + dze*dze);
        if (r_sat_earth <= EARTH_RADIUS) {
            printf("Crashed into Earth at t = %.1f s\n", time);
            break;
        }

        // RK4 step
        rk4_step(bodies, 3, DT);
        time += DT;

        if ((step+1) % 5000 == 0)
            printf("Step %d, time = %.1f days\n", step+1, time/86400.0);
    }

    fclose(fp);
    printf("\nSimulation finished.\n");
    if (loi_done) {
        printf("Lunar orbit insertion completed. Final orbit around Moon.\n");
    } else if (closest_approach < 1e30) {
        printf("Closest approach to Moon: %.1f km\n", closest_approach/1000.0);
    }
    printf("Results written to 'data/lunar_transfer.csv'\n");
    return 0;
}