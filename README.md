# Introduction

This is molecular dynamics sample program, which is able to simulate 1-atom quazy-ideal gas like Ar.

## Usage

To compile program, just run 

```
make
```

or simply run

```
./test.sh
```

and after compilation program runs automatically.

To compile program in debug mode, run

```
make debug
```
To delete executable file, run

```
make clean

```

## Modelling Argon

To start modelling, just type in UNIX terminal

```
./MD
```

You will see this output:
```
**********************
* Molecular Dynamics *
*  sample prorgram   *
*                    *
*       v 1.10       *
**********************
Enter number of steps:

```
So you need to enter number of steps, then you will have to enter number of particles and initial temperature:

```
**********************
* Molecular Dynamics *
*  sample prorgram   *
*                    *
*       v 1.10       *
**********************
Enter number of steps:
1000
Enter number of particles:
1000
Enter initial temperature in K:

```
This program is prepared for simulating Ar in normal conditions, for doing this, you have to enter  at least 1000 number of steps, number of particles should be equal to 1000 and temperature should be equal to 274:

```
**********************
* Molecular Dynamics *
*  sample prorgram   *
*                    *
*       v 1.10       *
**********************
Enter number of steps:
1000
Enter number of particles:
1000
Enter initial temperature in K:
274
```
After entering temperature, modelling will run:

```
**********************
* Molecular Dynamics *
*  sample prorgram   *
*                    *
*       v 1.10       *
**********************
Enter number of steps:
1000
Enter number of particles:
1000
Enter initial temperature in K:
30
Carrying out experiment...
Number of particles N: 1000
 Box size: Lx 3.35e-08 Ly 3.35e-08 Lz 3.35e-08
random generated T = 0.000400926
init T = 30
initial potential energy = 0
initial kinetic energy = 6.20828e-19
initial energy = 6.20828e-19
Density (m / V): 1.77415
Numerical density (N/V): 2.6599e+25
Theoretical pressure: 11008.9

```
You can obtain some data during run, but you should understand that system must be in termal equlibrium to make your data be correct. Let's wait while modelling is finished:

```
******************************** 
Modelling is finished, measuring final data
max_v: 220.977
step_v: 7.36589
max_prob_v N: 82
2nd max_prob_v N: 76
prob_v: 128.903
avg v: 130.664
Mass center impulse = 2.08537e-22
Kinetic Energy = 6.20828e-19
Full Energy = 6.20828e-19
Pressure = 10953.4 vs Theoretical pressure: 11008.9
Measuring current temperature 
Final T = 273.7917 vs Initial T = 274
******************************** 
```
As you can see form output above, the system density is 1.77415 in kg/m^3 and Temperature is about 274 K, so theory predicts, that the pressure of Ar is P = rho * k_B * T = 10^5 Pa. We can see, that our system's pressure is  10953.4, which is close enought to theoretical one, so, modelling went fine.

## Source code

This code is under MIT license

It is divided to 3 main files:
I. "main.cpp" - main file. It contains:
    1) namespace "ArgonSystem", which describes Ar gas. You can add your own namespace with different 1-molecule gas, for this you must specify mass of an atom, and Boltzman's const - they both will be used directly. Also you should implement calculating of potential, forces, and algoritms for calculating radii of atoms and velocities at every MD step. In "ArgonSystem" namespace we use Lennard-Jones potential (and its derivative as force) and Verlet algoritm for obtaining radii and velosities
    2)  main function with I/O and some variables. As you can see, it is possible to specify different dimension of system (or let it stay cubic), or change its timestep, termalization time (used for skipping initial part of trajectory for correct calculating of RDF and other things. For calculating RDF you should use external utils, for example, from Mdynamics:  http://www.fos.su.se), etc., but you should be careful - the Verlet algoritm can works badly just because it isn't suitable for this conditions (for example, timestep is too large).
II. "MDSys.cpp" - it contains MDSys class, which describes our system, and all methods we need to do a simulation. Class reads some data like number of atoms, mass of an atom, Boltzman's const, function for calcutating forces, etc - this allows you to use it with different gases, but there are essential limitations:
    * we assume that gas is 1-molecular (or this assumption is sutable for describing your system). No lipids and DNA are supported
    * we initialize system by placing atoms on a cubic lattice. Sorry, program can't read .XYZ or .xmol files
    * we output trajectory in .xmol format, to use different file format, you should implement writing yourself
III. "vector.h" - external library for 3d vectors math. You can write your own or use any librbay, which can calculate dot product, value of a vector, normalize it, etc, but notice that this methods are hardcoded in "MDSys.cpp".

code is briefly commented for better understanding how it works. You are free to use it for styding the Molecular Dynamics method

## Reading and useful links

http://www.fos.su.se - MDynamics main page

http://manual.gromacs.org/documentation/ - manual for GROMACS - another scientific MD tool

Understanding Molecular Simulation: From Algorithms to Applications, Book by Berend Smit and Daan Frenkel - useful book for studying molecular simulation

https://spiral.imperial.ac.uk/bitstream/10044/1/262/1/edm2006jcp.pdf - article about calculating pressure

