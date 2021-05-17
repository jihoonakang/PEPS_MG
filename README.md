## Poisson Equation Parallel Solver - Multigrid Method (PEPS_MG)

PEPS_MG provides a numerical solution of three-dimensional Poisson equation in Cartesian coordinate system with structured grids, employing parallel geometric multigrid methods. To handle the low granularity in coarse levels, the coarse-grid aggregation (CGA) procedure is implemented, which aggregates the grid variables at the certain level into a single core and proceeds the multigrid procedure until a sigle grid remains. A variation of CGA, named as coarse-grid partial semi-aggregation (CGPSA), is newly devised and implemented, where coarse grids in each coordinate dimension are selectively and independently step-by-step across multiple levels. With the CGPSA procedure, a communication structure becomes hierachycal and communication oeverheads are distributed into multiple levels and more computing cores can participate the grid processing on coarser levels.

## Features
- Three-dimensional Poisson equation with various boundary conditions (Needs to be implemented in matrix equations)
- Red-black Gauss-Siedel and conjugate gradient procedures as smoothers 
- Structured grids with unequal mesh sizes (Improvement required)
- MPI parallelization with three-dimensional domain decomposition and OpenMP parallelization for do-loops
- Three coarse grid solvers
  - (1) Parallel coarse grid solver - Obtain converged solution at the coarsest level.
  - (2) Coarse-grid aggregation - Aggregate grid variables into a single core at the aggregation level.
  - (3) Coarse-grid partial semi-aggregation - Aggregate grid variables step-by-step in hierachycal pattern.

## Authors
- Ji-Hoon Kang (jhkang@kisti.re.kr), Korea Institute of Science and Technology Information

## Cite
Ji-Hoon Kang, *"Scalable implementation of multigrid methods using partial semi-aggregation of coarse grids"*, under review.
