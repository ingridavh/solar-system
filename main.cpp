#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "solarsystem.h"
#include "solver.h"
#include <armadillo>
using namespace std;

int main(int numArguments, char **arguments)
{
    int numTimesteps = 2E08;
    if(numArguments >= 2) numTimesteps = atoi(arguments[1]);

    SolarSystem solarSystem;
    //-----------Celestial bodies---------------------------------

//    CelestialBody &earth = solarSystem.createCelestialBody( vec3(9.419288875250327E-01, 3.422743349115224E-01, -1.774653038679687E-04), vec3(-6.128263831462272E-03, 1.611761267097599E-02, 1.349643765318894E-07)*365, 3e-6 );
//    CelestialBody &jupiter = solarSystem.createCelestialBody(vec3(-5.429616996509673E+00, -4.392482185767863E-01, 1.232526518227290E-01), vec3(5.206515353907882E-04, -7.164038754622682E-03, 1.814205403145239E-05)*365., 0.001);
    CelestialBody &sun = solarSystem.createCelestialBody( vec3(0,0,0), vec3(0,0,0), 1.0 );
//    CelestialBody &mercury = solarSystem.createCelestialBody(vec3(-3.299771220462327E-01, 1.314596035268270E-01, 4.090462570196151E-02), vec3 (-1.589188578333726E-02, -2.505974994975800E-02, -5.904299543073971E-04)*365., 1.2E-07);
    //CelestialBody &mars = solarSystem.createCelestialBody(vec3 (1.145481237413500E+00, -7.759469057667274E-01, -4.451643705572472E-02), vec3 (8.418288645937397E-03, 1.276419831984888E-02, 6.072729084019031E-05)*365., 3.3E-07);
//    CelestialBody &venus = solarSystem.createCelestialBody(vec3 (1.619559463687653E-01, -7.064885525743501E-01, -1.903349841875333E-02), vec3 (1.960339728192615E-02, 4.341740156278193E-03, -1.071874641156530E-03)*365., 2.45E-06);
//    CelestialBody &saturn = solarSystem.createCelestialBody(vec3 (-2.277305999381926E+00, -9.772253310113138E+00, 2.605465118509321E-01), vec3 (5.126703359538530E-03, -1.283735869751270E-03, -1.816227490732958E-04)*365., 2.75E-4);
//    CelestialBody &uranus = solarSystem.createCelestialBody(vec3(1.846626259319022E+01, 7.554511790954641E+00, -2.111764845597773E-01), vec3(-1.517937917646755E-03, 3.456908604164539E-03, 3.251073912887524E-05)*365.,4.4E-05);
//    CelestialBody &neptune = solarSystem.createCelestialBody(vec3(2.825889439040351E+01, -9.928255088021324E+00, -4.468016759479966E-01), vec3(1.019309416815733E-03, 2.979917067133604E-03, -8.498151187073363E-05)*365., 5.15E-05);
//    CelestialBody &pluto = solarSystem.createCelestialBody(vec3(9.417687691553871E+00, -3.181872837345281E+01, 6.806467814675531E-01), vec3(3.073735631085241E-03, 2.467755912654963E-04, -9.027657565305481E-04)*365., 6.55E-09);


//    Give sun initial velocity so that sum (p) = 0
//    vec3 v_sun = jupiter.velocity*jupiter.mass*(-1) - earth.velocity*earth.mass;
//    CelestialBody &sun = solarSystem.createCelestialBody(vec3(3.569552387089993E-03, 3.395883113303513E-03, -1.598805663679657E-04), vec3(v_sun),1.0);

    //Perihelion
    CelestialBody &mercury = solarSystem.createCelestialBody(vec3(0.3075,0,0), vec3(0, 12.44,0), 1.2E-07);





    // To get a list (a reference, not copy) of all the bodies in the solar system, we use the .bodies() function
    vector<CelestialBody> &bodies = solarSystem.bodies();

    for(int i = 0; i<bodies.size(); i++) {
        CelestialBody &body = bodies[i]; // Reference to this body
        cout << "The position of this object is " << body.position << " with velocity " << body.velocity << endl;
    }

    double dt = 1E-7;
    Solver integrator(dt);

//    for(int timestep=0; timestep<numTimesteps; timestep++) {
//        integrator.Verlet(solarSystem);
//        solarSystem.writeToFile("positions_peri.txt");
//    }

//    cout << solarSystem.angularMomentum() << endl;
//    cout << "I just created my first solar system that has " << solarSystem.bodies().size() << " objects." << endl;


    // Set some helper variables before we start the time integration.
    double thetaPrevious 	= 0;	// The perihelion angle of the previous time step.
    double theta 		= 0;	// The perihelion angle of the current time step.

    double rPreviousPrevious 	= 0;	// Mercury-Sun-distance two times steps ago.
    double rPrevious   	 	= 0;	// Mercury-Sun-distance of the previous time step.
    double r 		 	= 0;	// Mercury-Sun-distance of the current time step.


//------------------Perihelion--------------------
    // This is the integration loop, in which you advance the solution (usually via a integrateOneStep()
    // function located in an integrator object, e.g. the Verlet class).
    //Write to file
    ofstream m_file;
    string filename = "perihelion.txt";
    if(m_file.good()) {
        m_file.open(filename, ofstream::out);
        if(!m_file.good()) {
            cout << "Error opening file " << filename << ". Aborting!" << endl;
            terminate();
        }
    }

    m_file << "Count of bodies " << solarSystem.numberOfBodies() << endl;
    for (int timeStep = 0; timeStep < numTimesteps; timeStep++) {
        integrator.Verlet(solarSystem);

        double x = bodies[1].position.x() - bodies[0].position.x();
        double y = bodies[1].position.y() - bodies[0].position.y();
        double thetaCurrent = atan2( y, x );

        double rCurrent = (bodies[1].position - bodies[0].position).length();

        if ( rCurrent > rPrevious && rPrevious < rPreviousPrevious ) {

            // If we are perihelion, print angle (in radians) to terminal.
            cout << thetaPrevious << endl;

            //Write results to file
            m_file << "Comment line that needs to be here. Blomst." << endl;
            m_file << "1 " << thetaPrevious << "\n";

        }

        // Update the helper variables (current, previous, previousPrevious).
        rPreviousPrevious 	= rPrevious;
        rPrevious		= rCurrent;
        thetaPrevious		= thetaCurrent;
    }
    return 0;
}

