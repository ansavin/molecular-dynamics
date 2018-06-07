#ifndef MDSYS_HPP
#define MDSYS_HPP

#include <fstream>

#include "vector.h"

#define PRINT_VEC(vec)                                                         \
    (cout << "(" << vec.getX() << ";" << vec.getY() << ";" << vec.getZ()       \
          << ")" << endl);
#define POTENTIAL_ENERGY_PERCENT 0.01
#define ENERGY_PERCENT 0.2

typedef double (*force_type)(double);
typedef double (*potential_type)(double);
typedef Vec3<double> vector;
typedef vector (*alg_R)(vector, vector, vector, double, double);
typedef vector (*alg_V)(vector, double);

class MDSys
{
  public:
    MDSys(force_type Force, alg_R Calc_R, alg_V Calc_V,
          potential_type Potential, double Mass, double K_B, double L_x,
          double L_y, double L_z, double PartNum, double tShift, double initT,
          double V0rand, double tTermal);
    ~MDSys();

    void initialize();
    void calcForces();
    void integrate();
    void outputData();

    std::ofstream out;      // output stream to .xmol file (use it for vmd)
    std::ofstream term_out; // output stream to .001 file (use it for RDF)

  private:
    void measureValues();
    double calcMaxProbabV() const;
    double calcPressure();
    void printToFile(int i);
    void newFrame();
    void substractSysVelocity();
    vector NIM(vector r1, vector r2, double Dx, double Dy, double Dz) const;
    double shiftNIM(double coord, double l) const;
    void PBC(vector &dist, double Dx, double Dy, double Dz) const;
    double correctCoord(double coord, double leftBoundary,
                        double rightBoundary) const;

    int N;                         // number of the particles
    double Lx = 0, Ly = 0, Lz = 0; // box size in meters
    double dt;                     // MD time step in seconds
    double T0;                     // Initial temperature in K
    bool v0_rand;                  // do particles have initial velocities?
    bool isInitialised = false;    // is system initialised now?
    vector *r;                     // all particles coordinates
    vector *dr;    // distance from old to current coords for all particles
    vector *v;     // velocities of all particles
    vector *f;     // forces acting on a particle for every particle
    double m = 0;  // mass of particle
    double k_B;    // Boltzman`s const
    double energy; // system`s full energy
    double V;      // system`s volume
    double rho;    // system`s density
    double initial_energy; // initial energy of the system
    double
        L_FREE_MOTION; // const for checking that r(t+dt) - r(t) isn`t too large
    double potential_energy; // potenial energy of the system
    double kinetic_energy;   // kinetic energy of the system
    double averV;            // average velocity of system`s particle
    double sqrAverV;         // squared average velocity of system`s particle
    double T;                // system`s current temerature
    double pressure;         // system`s pressure
    vector momentum;         // system`s momentum
    force_type force;        // force acting between patricles
    alg_R calc_R;            // algoritm for calculating radii of particles
    alg_V calc_V;            // algoritm for calculating velocities of particles
    potential_type potential;  // potential acting between patricles
    double time;               // current modelling time
    double termalization_time; // time for system to termalise in sec
};

#endif
