// uncomment to disable assert()
//#define NDEBUG
#include <cassert>

#include <ctime>

#include <fstream>
#include <iomanip>
#include <memory.h>
#include <omp.h>
#include <sstream>

#include "MDSys.h"
#include "vector.h"

MDSys::MDSys(force_type Force, alg_R Calc_R, alg_V Calc_V,
             potential_type Potential, double Mass, double K_B, double L_x,
             double L_y, double L_z, double PartNum, double tShift,
             double initT, double V0rand, double tTermal)

    : N(PartNum), Lx(L_x), Ly(L_y), Lz(L_z), dt(tShift), T0(initT),
      v0_rand(V0rand), m(Mass), k_B(K_B), force(Force), calc_R(Calc_R),
      calc_V(Calc_V), potential(Potential), termalization_time(tTermal)
{
    dr = new vector[N]; // distance from old to current coords
    v = new vector[N];  // velocities
    f = new vector[N];  // forces acting on a particle
    r = new vector[N];  // current coordinates

    time = 0;           // time of current MD step
    energy = 0;         // full energy
    kinetic_energy = 0; // kinetic energy
    initial_energy = 0; // initial full energy

    out.open("output.xmol");     // output file
    term_out.open("output.001"); // output file for RDF

    // we use it to check that r(t+dt) - r(t) isn`t too large
    L_FREE_MOTION = (pow(Lx * Ly * Lz, 1.0 / 3) / (2 * N));
    V = Lx * Ly * Lz;
}

MDSys::~MDSys()
{
    delete[] dr;
    delete[] v;
    delete[] f;
    delete[] r;

    out.close();
    term_out.close();
}

/*
 * function for putting particles on a lattice
 * and giving them initial speeds and coordinates
 * in previous time (for correct work of Verlet
 * algoritm)
 */
void MDSys::initialize()
{
    srand(std::time(NULL));

    int K = ceil(pow(N, 1.0 / 3));
    double dLx = Lx / K;
    double dLy = Ly / K;
    double dLz = Lz / K;

    newFrame();
    std::cout << "Number of particles N: " << N << std::endl;
    std::cout << " Box size: Lx " << Lx << " Ly " << Ly << " Lz " << Lz
              << std::endl;

    if (v0_rand)
    {
        // ok, we want particles to have initial random velocities
        for (int i = 0; i < N; i++)
        {
            v[i] = vector((rand() % 10000 - 5000) / 10000.0,
                          (rand() % 10000 - 5000) / 10000.0,
                          (rand() % 10000 - 5000) / 10000.0);
        }

        // calc temperature and so on...
        measureValues();

        std::cout << "random generated T = " << T << std::endl;

        // fix temperature to needed
        for (int i = 0; i < N; i++)
        {
            v[i] = vector(v[i].getX() * pow(T0 / T, 0.5),
                          v[i].getY() * pow(T0 / T, 0.5),
                          v[i].getZ() * pow(T0 / T, 0.5));
        }

        measureValues();

        std::cout << "init T = " << T << std::endl;

        // ok, we want system at all to stay at rest
        substractSysVelocity();
    }
    else
    {
        for (int i = 0; i < N; i++)
        {
            // ok, we want particles to stay in rest
            v[i] = vector(0, 0, 0);
        }
    }

    for (int counter = 0, i = 0; counter < N && i < K; i++)
    {
        for (int j = 0; counter < N && j < K; j++)
        {
            for (int k = 0; counter < N && k < K; k++)
            {
                // put particle into lattice
                r[counter].set((i + 1.0 / 2) * dLx, (j + 1.0 / 2) * dLy,
                               (k + 1.0 / 2) * dLz);
                // init dr for Verle algo step
                dr[counter].set(v[counter].getX() * 2 * dt,
                                v[counter].getY() * 2 * dt,
                                v[counter].getZ() * 2 * dt);

                // check that particle move not too fast
                assert(dr[counter].abs() < L_FREE_MOTION);

                printToFile(counter);

                counter++;
            }
        }
    }

    isInitialised = true;

    measureValues();
    initial_energy = energy;

    std::cout << "initial potential energy = " << potential_energy << std::endl;
    std::cout << "initial kinetic energy = " << kinetic_energy << std::endl;
    std::cout << "initial energy = " << energy << std::endl;

    // check ideality
    assert(!v0_rand || (fabs(potential_energy / initial_energy) <
                        POTENTIAL_ENERGY_PERCENT));

    rho = m * N / V;
    std::cout << "Density (m / V): " << rho << std::endl;
    std::cout << "Numerical density (N/V): " << rho / m << std::endl;
    std::cout << "Theoretical pressure: " << rho * k_B * T / m << std::endl;
}

/*
 * Method for calculating forces acting on all particles
 */
void MDSys::calcForces()
{
    double ff;
    double _rij;

#pragma omp parallel for private(ff, _rij)
    for (int i = 0; i < N; i++)
    {
        f[i] = vector(0, 0, 0);
        for (int j = 0; j < N; j++)
        {
            if (i != j)
            {
                vector rij = NIM(r[i], r[j], Lx, Ly, Lz); // r[i] - r[j]
                _rij = rij.abs();                         // |r[i] - r[j]|
                ff = force(_rij);
                vector dr = rij;
                dr.normalise();
                f[i] += dr * ff;
            }
        }
    }
}

/*
 * Method integrate differential equations
 * and make one step modelling algo
 */
void MDSys::integrate()
{
    vector rBuf;

    time += dt;

    newFrame();

#pragma omp parallel for private(rBuf)
    for (int i = 0; i < N; i++)
    {
        rBuf = r[i];

        // calculating coordinates using external algorithm
        r[i] = calc_R(r[i], dr[i], f[i], m, dt);
        dr[i] = r[i] - rBuf;

        // check that particle move not too fast
        assert(dr[i].abs() < L_FREE_MOTION);

        // calculate velocities using external algo
        v[i] = calc_V(dr[i], dt);

        // use PBC method
        PBC(r[i], Lx, Ly, Lz);
    }

    // dump data to file
    for (int i = 0; i < N; i++)
    {
        printToFile(i);
    }
}

/*
 * function for measuring and output some variables
 */
void MDSys::outputData()
{

    measureValues();
    double max_prob_v = calcMaxProbabV();
    std::cout << "prob_v: " << max_prob_v << std::endl;
    std::cout << "avg v: " << averV << std::endl;

    std::cout << "Mass center impulse = " << momentum.abs() << std::endl;
    // ok, we want system at all to stay in rest during modelling
    substractSysVelocity();

    std::cout << "Kinetic Energy = " << kinetic_energy << std::endl;
    std::cout << "Full Energy = " << energy << std::endl;
    std::cout << "Pressure = " << pressure << " vs "
              << "Theoretical pressure: " << rho * k_B * T / m << std::endl;
    std::cout << "Measuring current temperature " << std::endl;
    std::cout << "Final T = " << T << " vs "
              << "Initial T = " << T0 << std::endl;
}

/*
 * Writes some simulation stage data to out stream
 */
void MDSys::newFrame()
{
    out << N << std::endl;
    out << "***** time = " << time << " ***** " << std::endl;
    if (time > termalization_time)
    {
        term_out << N << std::endl;
        term_out << "***** time = " << time << " ***** " << std::endl;
    }
}

/*
 * Method for substracting system`s speed from particle velocities
 * to make system at all stay in rest
 */

void MDSys::substractSysVelocity()
{
#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        v[i] -= momentum / (N * m);
    }
}

/*
 * print to file data of particle #i
 */
void MDSys::printToFile(int i)
{
    out << "Ar " << std::setw(10) << r[i].getX() << " " << std::setw(10)
        << r[i].getY() << " " << std::setw(10) << r[i].getZ() << " "
        << std::setw(10) << v[i].getX() << " " << std::setw(10) << v[i].getY()
        << " " << std::setw(10) << v[i].getZ() << std::endl;
    if (time > termalization_time)
    {
        term_out << "Ar " << std::setw(10) << r[i].getX() << " "
                 << std::setw(10) << r[i].getY() << " " << std::setw(10)
                 << r[i].getZ() << " " << std::setw(10) << v[i].getX() << " "
                 << std::setw(10) << v[i].getY() << " " << std::setw(10)
                 << v[i].getZ() << std::endl;
    }
}

/*
 * Nearest image method implementation
 */
vector MDSys::NIM(vector r1, vector r2, double Dx, double Dy, double Dz) const
{
    double y = -(r1.getY() - r2.getY());
    double x = -(r1.getX() - r2.getX());
    double z = -(r1.getZ() - r2.getZ());

    vector dist = vector(0, 0, 0);

    dist.setX(shiftNIM(x, Dx));
    dist.setY(shiftNIM(y, Dy));
    dist.setZ(shiftNIM(z, Dz));

    // check that after distance normalization we have it less than rectangle
    // diagonal
    assert(dist.abs() < vector(Dx, Dy, Dz).abs());

    return dist;
}

/*
 * Method that shift particle coordinates for nearest image method
 * implementation
 */
double MDSys::shiftNIM(double coord, double l) const
{
    if (coord >= l / 2.0)
    {
        coord = l - coord;
    }
    else if (coord <= -l / 2.0)
    {
        coord = coord + l;
    }

    return coord;
}

/*
 * Periodic boundary conditions implementation
 */
void MDSys::PBC(vector &dist, double Dx, double Dy, double Dz) const
{
    double aDx = 0;
    double aDy = 0;
    double aDz = 0;

    dist.setX(correctCoord(dist.getX(), aDx, Dx));
    dist.setY(correctCoord(dist.getY(), aDy, Dy));
    dist.setZ(correctCoord(dist.getZ(), aDz, Dz));
}

/*
 * Periodic boundary conditions coordinates correction implementation
 */
double MDSys::correctCoord(double coord, double leftBoundary,
                           double rightBoundary) const
{
    double l = rightBoundary - leftBoundary;
#ifndef NDEBUG
    std::stringstream dbg;
    if (coord >= rightBoundary || coord < leftBoundary)
    {
        dbg << "Function PBC starts,  " << coord << " !in "
            << "(" << leftBoundary << ";" << rightBoundary << ")" << std::endl;
    }
#endif

    double d;

    if (coord >= rightBoundary)
    {
        d = (coord - leftBoundary);
#ifndef NDEBUG
        dbg << "Crossing right boundary:: (d / l) : floor(d / l) : (l *  "
               "floor ( d / l ))"
            << (d / l) << " : " << floor(d / l) << " ; " << (l * floor(d / l))
            << std::endl;
#endif
        coord = (coord - l * floor(d / l));
    }
    else if (coord < leftBoundary)
    {
        d = (leftBoundary - coord);
#ifndef NDEBUG
        dbg << "Crossing left boundary:: (d / l) : floor(d / l) : (l *  "
               "floor ( d / l ))"
            << (d / l) << " : " << floor(d / l) << " ; "
            << (d - l * floor(d / l)) << std::endl;
#endif
        coord = (rightBoundary - l * ((d / l) - floor(d / l)));
    }

#ifndef NDEBUG
    if (coord > rightBoundary || coord < leftBoundary)
    {
        // dbg >> std::cout;
        std::cout << dbg.rdbuf() << "Function falure " << coord << " !in "
                  << "(" << leftBoundary << ";" << rightBoundary << ")"
                  << std::endl;
        exit(-1);
    }
#endif

    return coord;
}

/*
 * calculating some physical values like average v
 */

void MDSys::measureValues()
{

    double avgV = 0;
    double sqrAvgV = 0;
    double kinEnergy = 0;
    double potEnergy = 0;
    double E = 0;
    const double cor_ratio = sqrt(3 * 3.1416 / 8);
    double ratio = 0;
    momentum = vector(0, 0, 0);
    pressure = 0;

#pragma omp parallel for reduction(+ : avgV, sqrAvgV, kinEnergy, potEnergy)
    for (int i = 0; i < N; i++)
    {
        kinEnergy += m * v[i].abs2() / 2;
        momentum += v[i] * m;
        if (isInitialised)
        {
            avgV += v[i].abs();
            sqrAvgV += v[i].abs2();
            for (int j = 0; j < N; j++)
                if (i != j)
                {
                    vector dr = NIM(r[i], r[j], Lx, Ly, Lz); // r[i] - r[j];
                    double rij = dr.abs();
                    potEnergy += potential(rij);
                }
        }
        pressure += m * v[i].abs2() + vector::dotProduct(r[i], f[i]);
    }
    avgV /= N;
    sqrAvgV = sqrt(sqrAvgV / N);
    potEnergy /= 2;
    T = kinEnergy / (3.0 * N * k_B / 2);
    E = potEnergy + kinEnergy;
    pressure /= 3 * V;
    /*
    asserts to ensure that modelling is correct
    */
    if (isInitialised)
    {
        ratio = sqrAvgV / avgV;
        assert(fabs(ratio - cor_ratio) / cor_ratio < 0.05);

        // energy conservation check
        // assert(fabs(E - initial_energy) < ENERGY_PERCENT *
        // fabs(initial_energy));
    }
    averV = avgV;
    sqrAverV = sqrAvgV;
    kinetic_energy = kinEnergy;
    potential_energy = potEnergy;
    energy = E;
}

/*
 * Debug purpose method for check modelling correctness by
 * analyzing max probable velocity
 */
double MDSys::calcMaxProbabV() const
{
    const int intervals = 30;
    double max_v = 0;

    int *v_list = new int[intervals];
    memset(v_list, 0, sizeof(int) * intervals);
    double tmp_v;

    for (int i = 0; i < N; i++)
    {
        tmp_v = v[i].abs();
        if (max_v < tmp_v)
            max_v = tmp_v;
    }

    double delta = max_v / intervals;

    int max_group = 0;
    int max_group_prob = 0;
    int _2nd_max_group = 0;
    int _2nd_max_group_prob = 0;
    for (int i = 0; i < N; i++)
    {
        tmp_v = v[i].abs();
        int group_no = static_cast<int>(floor(tmp_v / delta));
        group_no = (group_no < intervals) ? group_no : (group_no - 1);

        v_list[group_no]++;

        if (v_list[group_no] > max_group_prob)
        {
            max_group = group_no;
            max_group_prob = v_list[group_no];
        }
        else if (v_list[group_no] > _2nd_max_group_prob)
        {
            _2nd_max_group = group_no;
            _2nd_max_group_prob = v_list[group_no];
        }
    }

    std::cout << "max_v: " << max_v << std::endl;
    std::cout << "step_v: " << delta << std::endl;
    std::cout << "max_prob_v N: " << v_list[max_group] << std::endl;
    std::cout << "2nd max_prob_v N: " << v_list[_2nd_max_group] << std::endl;

    delete[] v_list;

    return (max_group + 0.5) * delta;
}
