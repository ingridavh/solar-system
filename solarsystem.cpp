#include "solarsystem.h"
#include <iostream>
#include <cmath>
using namespace std;

SolarSystem::SolarSystem() :
    m_kineticEnergy(0),
    m_potentialEnergy(0),
    m_G(4*M_PI*M_PI),
    m_c(63072.0)
{

}

CelestialBody& SolarSystem::createCelestialBody(vec3 position, vec3 velocity, double mass) {
    m_bodies.push_back( CelestialBody(position, velocity, mass) );
    return m_bodies.back(); // Return reference to the newest added celstial body
}

void SolarSystem::calculateForcesAndEnergy()
{
    m_kineticEnergy = 0;
    m_potentialEnergy = 0;
    m_angularMomentum.zeros();

    for(CelestialBody &body : m_bodies) {
        // Reset forces on all bodies
        body.force.zeros();
    }

    for(int i=0; i<numberOfBodies(); i++) {
        CelestialBody &body1 = m_bodies[i];
        for(int j=i+1; j<numberOfBodies(); j++) {
            CelestialBody &body2 = m_bodies[j];
            vec3 deltaRVector = body1.position - body2.position;
            vec3 deltaVvector = body1.velocity - body2.velocity;
            double dr = deltaRVector.length();
            //Calculate potential energy
            m_potentialEnergy -= m_G*body1.mass*body2.mass/dr;

            //Calculate the force on both bodies
            //body1.force -= (m_G*body1.mass*body2.mass)/(dr*dr*dr)*deltaRVector;
            //body2.force += (m_G*body1.mass*body2.mass)/(dr*dr*dr)*deltaRVector;

            //Calculate relativistic correction to force
            double l = (deltaRVector.cross(deltaVvector)).length();
            body1.force -= (m_G*body1.mass*body2.mass)/(dr*dr*dr)*(1 + 3*l*l/(dr*dr*m_c*m_c))*deltaRVector;
            body2.force += (m_G*body1.mass*body2.mass)/(dr*dr*dr)*(1 + 3*l*l/(dr*dr*m_c*m_c))*deltaRVector;


            //Calculate potential energy for each body
            body1.potential_energy -= m_G*body1.mass*body2.mass/dr;
            body2.potential_energy -= m_G*body1.mass*body2.mass/dr;


        }

        body1.kinetic_energy = 0.5*body1.mass*body1.velocity.lengthSquared();
        m_kineticEnergy += body1.kinetic_energy;

        //Why does L mess everything up?
        m_angularMomentum = body1.mass*body1.position.cross(body1.velocity);

        //Perihelion force


    }
}

int SolarSystem::numberOfBodies() const
{
    return m_bodies.size();
}

double SolarSystem::totalEnergy() const
{
    return m_kineticEnergy + m_potentialEnergy;
}

double SolarSystem::potentialEnergy() const
{
    return m_potentialEnergy;
}

double SolarSystem::kineticEnergy() const
{
    return m_kineticEnergy;
}

void SolarSystem::writeToFile(string filename)
{
    if(!m_file.good()) {
        m_file.open(filename.c_str(), ofstream::out);
        if(!m_file.good()) {
            cout << "Error opening file " << filename << ". Aborting!" << endl;
            terminate();
        }
    }

    m_file << "Count of bodies " << numberOfBodies() << endl;
    m_file << "Comment line that needs to be here. Blomst." << endl;
    for(CelestialBody &body : m_bodies) {
        m_file << "1 " << body.position.x() << " " << body.position.y() << " " << body.position.z() << "\n";
    }
}

vec3 SolarSystem::angularMomentum() const
{
    return m_angularMomentum;
}

std::vector<CelestialBody> &SolarSystem::bodies()
{
    return m_bodies;
}

double SolarSystem::G() const
{
    return m_G;
}

void SolarSystem::setG(double G)
{
    m_G = G;
}



















