# ExtendedLagrangian

Extended Lagrangian Schemes for Solving Charge Equilibrium in ReaxFF

The SC-XLMD scheme involves the modification of 7 files in LAMMPS, which are listed in `sc-xlmd/`. The original file, from the Aug 2019 release of LAMMPS, is placed in `original/` for comparison.

## Diff code

You can use `diff` to compare the ones in `sc-xlmd/` and `original/`. For example, use

```
diff sc-xlmd/atom.h original/atom.h
```

## Build LAMMPS

The function `Langevin()` uses C++ 2011 standard library to generate normal random numbers, so you need to build with standard `-std=c++11`, or rewrite it in other ways to avoid this.

```
cd build
cmake -D PKG_USER-REAXC=1 -D CMAKE_CXX_COMPILER=mpiicpc -D CMAKE_CXX_FLAGS='-std=c++11' ../lammps/cmake > cmake.log
make -j8 > make.log
```

## Run an example

After you build the LAMMPS binary, you can use

```
mpiexec -n 8 build/lmp -i assets/water.in -l assets/water.log -sc none -nc
```
