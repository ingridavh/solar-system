#ifndef SOLVER_H
#define SOLVER_H


class Solver
{
public:
    double m_dt;
    Solver(double dt);
    void Euler(class SolarSystem &system);
    void Verlet(class SolarSystem &system);
//    void integrateOneStep(class SolarSystem &system);
};

#endif // SOLVER_H
