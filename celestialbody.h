#ifndef CELESTIALBODY_H
#define CELESTIALBODY_H

#include "vec3.h"

class CelestialBody
{
public:
    vec3 position;
    vec3 velocity;
    vec3 force;
    double mass;
    double kinetic_energy;
    double potential_energy;
    double total_energy;

    CelestialBody(vec3 position, vec3 velocity, double mass);
    CelestialBody(double x, double y, double z, double vx, double vy, double vz, double mass);
    void resetForce();
    void resetEnergy();
};

#endif // CELESTIALBODY_H
