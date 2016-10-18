#include "solver.h"
#include "solarsystem.h"



Solver::Solver(double dt) :
    m_dt(dt)
{

}

void Solver::Euler(SolarSystem &system)
{
    //Calculate force on bodies in the system
    system.calculateForcesAndEnergy();


    //Integrate one step with the Euler method
    for(CelestialBody &body : system.bodies()) {
        body.position += body.velocity*m_dt;
        body.velocity += body.force / body.mass * m_dt;
    }
}


void Solver::Verlet(SolarSystem &system)
{
    //Calculate force on bodies in the system
    system.calculateForcesAndEnergy();

    //First intermediate step with a_i
    for(CelestialBody &body : system.bodies()) {
        body.velocity += 0.5*body.force/body.mass*m_dt;
        body.position += body.velocity*m_dt;
    }

    //Calculate new a_i+1
    system.calculateForcesAndEnergy();

    //Second intermediate step using a_i+1
    for(CelestialBody &body : system.bodies()) {
        body.velocity += 0.5*body.force / body.mass*m_dt;
    }
}
