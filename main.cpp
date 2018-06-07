#include "MDSys.h"
#include "vector.h"

#include <ctime>
#include <fenv.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <memory.h>
#include <omp.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

namespace ArgonSystem // namespace for modelling Argon in SI unit system
{
const double eps = 1.654e-21; // eps parameter of Ar for LD potential in Jowels
const double sigma =
    3.405e-10; // sigma parameter of Ar for LD potential in meters
const double rcut =
    2.5 * sigma; // standart rcut - if |r_i - r_j| > rcut then f = 0
const double rmin =
    0.00001 * sigma;         // if |r_i - r_j| < rmin then we count f using rmin
const double k_B = 1.38e-23; // the Boltsmann`s constant
const double mass = 6.67e-26; // the mass of an simgle Ar atom
/*
Lennard-Jones force
*/
double forceLD(double r) //+
{
    if (r > rcut)
        return 0;

    if (r < rmin)
        return forceLD(rmin);

    double x = sigma / r;
    return -48 * eps / sigma * (pow(x, 13) - 0.5 * pow(x, 7));
}

/*
Lennard-Jones potential
*/
double potentialLD(double r) //+
{
    if (r > rcut)
        return 0;

    if (r < rmin)
        return potentialLD(rmin);

    double x = sigma / r;
    return 4 * eps * (pow(x, 12) - pow(x, 6));
}

/*
Verlet algoritm:
*/
vector VerleR(vector r, vector dr, vector f, double m, double dt) //+
{

    return r + (dr + (f / (2 * m)) * dt * dt);
}

vector VerleV(vector dr, double dt) //+
{
    return dr / (2 * dt);
}

} // end of namespace

/*
 * Main function
 */
int main()
{
    // http://en.cppreference.com/w/cpp/numeric/math/math_errhandling
    feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

    omp_set_num_threads(8); // set the number of threads

    cout << "**********************" << endl;
    cout << "* Molecular Dynamics *" << endl;
    cout << "*  sample prorgram   *" << endl;
    cout << "*                    *" << endl;
    cout << "*       v 1.10       *" << endl;
    cout << "**********************" << endl;

    int numberOfSteps, particlesNum;
    double initT;

    // entering some variables
    cout << "Enter number of steps:" << endl;
    cin >> numberOfSteps;
    cout << "Enter number of particles:" << endl;
    cin >> particlesNum;
    cout << "Enter initial temperature in K:" << endl;
    cin >> initT;

    /*
     * This variables were well chosen for correct work
     * of the program. You can change it, but in this case
     * some problems may appear - just because the algoritms
     * used here may not be suitable for new conditions. This
     * is why this numbers are hardcoded
     */
    const double L_x = 3.35e-8, // box size in meters
        L_y = L_x, L_z = L_x;
    const double tShift = 1e-15; // MD time step in seconds
    bool V0rand = true;          // we need random initial speeds
    double tTermal = 1e-14;      // termalization time

    MDSys sample(ArgonSystem::forceLD, ArgonSystem::VerleR, ArgonSystem::VerleV,
                 ArgonSystem::potentialLD, ArgonSystem::mass, ArgonSystem::k_B,
                 L_x, L_y, L_z, particlesNum, tShift, initT, V0rand, tTermal);
    cout << "Carrying out experiment..." << endl;
    sample.initialize(); // initialise system

    // we want to print data 10 times during modelling
    int timeForKinECount = 0.1 * numberOfSteps;
    int counterForKinE = 0;

    double start = omp_get_wtime(); // start calculating time

    for (int i = 0; i < numberOfSteps; i++)
    {
        sample.calcForces(); // calc forces acting in our system
        sample.integrate();  // integrate equation of motion
        if (counterForKinE == timeForKinECount)
        {
            cout << "******************************** " << endl;
            cout << "MD step # " << i
                 << " stage modelling time: " << (omp_get_wtime() - start)
                 << endl;
            sample.outputData(); // measure data and averages

            counterForKinE = 0;
            start = omp_get_wtime();
        }
        else
        {
            counterForKinE++;
        }
    }
    // measure final data
    cout << "******************************** " << endl;
    cout << "Modelling is finished, measuring final data" << endl;
    sample.outputData();
    cout << "******************************** " << endl;

    return 0;
}
