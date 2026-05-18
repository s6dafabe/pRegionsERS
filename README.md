# pRegionsERS
Code for the paper "Strong ILP Formulations for the p-Regions Problem"

## Dependencies
The project has the following dependencies:
  - Boost 1.89
  - Gurobi 12

## Build instructions
```
  mkdir build
  cd build
  cmake ..
  make
```
## Usage
Use like this:
```
  ./aggregation <path/to/instance> <k> <Flow:Tree:ERS:ERSTree> <timelimit>
```
where k is the number of output regions and one of Flow|Tree|ERS|ERSTree determines the tested model

Example:
```
  ./aggregation instances/bulgaria.adj 4 ERSTree 3600
```

