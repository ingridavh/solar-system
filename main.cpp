#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "solarsystem.h"
#include "solver.h"
#include <armadillo>
using namespace std;

int main(int numArguments, char **arguments)
{
    int numTimesteps = 1500;
    if(numArguments >= 2) numTimesteps = atoi(arguments[1]);

    SolarSystem solarSystem;
    // We create new bodies like this. Note that the createCelestialBody function returns a reference to the newly created body
    // This can then be used to modify properties or print properties of the body if desired
    // Use with: solarSystem.createCelestialBody( position, velocity, mass );

    CelestialBody &sun = solarSystem.createCelestialBody( vec3(0,0,0), vec3(0,0,0), 1.0 );

    // We don't need to store the reference, but just call the function without a left hand side
    solarSystem.createCelestialBody( vec3(9.419288875250327E-01, 3.422743349115224E-01, -1.774653038679687E-04), vec3(-6.128263831462272E-03, 1.611761267097599E-02, 1.349643765318894E-07)*365., 3e-6 );

    //CelestialBody &jupiter = solarSystem.createCelestialBody(vec3(-5.429616996509673E+00,-4.392482185767863E-01,1.232526518227290E-01), vec3(5.206515353907882E-04,-7.164038754622682E-03,1.814205403145239E-05)*365., 0.001);


    // To get a list (a reference, not copy) of all the bodies in the solar system, we use the .bodies() function
    vector<CelestialBody> &bodies = solarSystem.bodies();

    for(int i = 0; i<bodies.size(); i++) {
        CelestialBody &body = bodies[i]; // Reference to this body
        cout << "The position of this object is " << body.position << " with velocity " << body.velocity << endl;
    }

    double dt = 0.001;
    Solver integrator(dt);

    //Create matrix to store energies [kinetic, potential, total]
    arma::mat dE_mat = arma::zeros<arma::mat>(numTimesteps,3);


    for(int timestep=0; timestep<numTimesteps; timestep++) {
        integrator.Verlet(solarSystem);
        //solarSystem.writeToFile("positions_jup_1000.txt");
        cout << solarSystem.angularMomentum() << endl;

        //Save energies in matrix
        dE_mat(timestep,0) = solarSystem.totalEnergy();
        dE_mat(timestep,1) = solarSystem.kineticEnergy();
        dE_mat(timestep,2) = solarSystem.potentialEnergy();
    }

    cout << "The maximum error for total energy conservation is " << abs ((dE_mat.col(0).max() - dE_mat.col(0).min())/dE_mat.col(0).max()) << endl;
    cout << "The maximum error for kinetic energy conservation is "<< abs ((dE_mat.col(1).max() - dE_mat.col(1).min())/dE_mat.col(1).max())<< endl;
    cout << "The maximum error for potenial energy conservation is "<< abs ((dE_mat.col(2).max() - dE_mat.col(2).min())/dE_mat.col(2).max()) << endl;
    cout << "I just created my first solar system that has " << solarSystem.bodies().size() << " objects." << endl;
    return 0;
}

