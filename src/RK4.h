#ifndef RK4_H
#define RK4_H

#include <math.h>

#ifndef G
#define G 6.67430e-11

/* Masses (kg) */
#define EARTH_MASS   5.972e24
#define MOON_MASS    7.34767309e22
#define SATELLITE_MASS 1000.0

/* Radii (m) */
#define EARTH_RADIUS 6371000.0
#define MOON_RADIUS  1737100.0
#endif

typedef struct {
    double x, y, z;
    double vx, vy, vz;
    double mass;
} Body;

typedef struct {
    double dx, dy, dz;
    double dvx, dvy, dvz;
} Derivative;


void compute_acceleration(Body *bi, Body *bj, double *ax, double *ay, double *az) {
    double dx = bj->x - bi->x;
    double dy = bj->y - bi->y;
    double dz = bj->z - bi->z;
    double distSq = dx*dx + dy*dy + dz*dz + 1e-10;
    double dist = sqrt(distSq);
    double force = G * bj->mass / (distSq * dist);
    *ax += force * dx;
    *ay += force * dy;
    *az += force * dz;
}

void compute_derivatives(Body bodies[], Derivative derivs[], int n) {
    for (int i = 0; i < n; i++) {
        derivs[i].dx = bodies[i].vx;
        derivs[i].dy = bodies[i].vy;
        derivs[i].dz = bodies[i].vz;
        double ax = 0, ay = 0, az = 0;
        for (int j = 0; j < n; j++)
            if (i != j)
                compute_acceleration(&bodies[i], &bodies[j], &ax, &ay, &az);
        derivs[i].dvx = ax;
        derivs[i].dvy = ay;
        derivs[i].dvz = az;
    }
}

void rk4_step(Body bodies[], int n, double dt) {
    Derivative k1[n], k2[n], k3[n], k4[n];
    Body temp[n];
    compute_derivatives(bodies, k1, n);
    for (int i = 0; i < n; i++) {
        temp[i] = bodies[i];
        temp[i].x += 0.5 * dt * k1[i].dx;
        temp[i].y += 0.5 * dt * k1[i].dy;
        temp[i].z += 0.5 * dt * k1[i].dz;
        temp[i].vx += 0.5 * dt * k1[i].dvx;
        temp[i].vy += 0.5 * dt * k1[i].dvy;
        temp[i].vz += 0.5 * dt * k1[i].dvz;
    }
    compute_derivatives(temp, k2, n);
    for (int i = 0; i < n; i++) {
        temp[i] = bodies[i];
        temp[i].x += 0.5 * dt * k2[i].dx;
        temp[i].y += 0.5 * dt * k2[i].dy;
        temp[i].z += 0.5 * dt * k2[i].dz;
        temp[i].vx += 0.5 * dt * k2[i].dvx;
        temp[i].vy += 0.5 * dt * k2[i].dvy;
        temp[i].vz += 0.5 * dt * k2[i].dvz;
    }
    compute_derivatives(temp, k3, n);
    for (int i = 0; i < n; i++) {
        temp[i] = bodies[i];
        temp[i].x += dt * k3[i].dx;
        temp[i].y += dt * k3[i].dy;
        temp[i].z += dt * k3[i].dz;
        temp[i].vx += dt * k3[i].dvx;
        temp[i].vy += dt * k3[i].dvy;
        temp[i].vz += dt * k3[i].dvz;
    }
    compute_derivatives(temp, k4, n);
    for (int i = 0; i < n; i++) {
        bodies[i].x += dt / 6.0 * (k1[i].dx + 2*k2[i].dx + 2*k3[i].dx + k4[i].dx);
        bodies[i].y += dt / 6.0 * (k1[i].dy + 2*k2[i].dy + 2*k3[i].dy + k4[i].dy);
        bodies[i].z += dt / 6.0 * (k1[i].dz + 2*k2[i].dz + 2*k3[i].dz + k4[i].dz);
        bodies[i].vx += dt / 6.0 * (k1[i].dvx + 2*k2[i].dvx + 2*k3[i].dvx + k4[i].dvx);
        bodies[i].vy += dt / 6.0 * (k1[i].dvy + 2*k2[i].dvy + 2*k3[i].dvy + k4[i].dvy);
        bodies[i].vz += dt / 6.0 * (k1[i].dvz + 2*k2[i].dvz + 2*k3[i].dvz + k4[i].dvz);
    }
}

#endif