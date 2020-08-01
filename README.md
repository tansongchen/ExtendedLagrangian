# ExtendedLagrangian

Extended Lagrangian Schemes for Solving Charge Equilibrium in ReaxFF

The SC-XLMD scheme involves the modification of 7 files in LAMMPS, which are listed in `sc-xlmd/`. The original file, from the Aug 2019 release of LAMMPS, is placed in `original/` for comparison.

You can use `diff` to compare the ones in `sc-xlmd/` and `original/`. For example, use

```
diff sc-xlmd/atom.h original/atom.h
```