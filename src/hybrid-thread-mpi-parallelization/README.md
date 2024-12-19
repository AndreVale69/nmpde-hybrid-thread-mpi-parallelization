# Hybrid thread-MPI parallelization

## Problem Description

Distributed memory parallelization is convenient when the different computing cores
are connected to different memory units. Conversely, when multiple cores have access
to the same memory, shared memory parallelization (multithreading) can be more
convenient, to minimize memory redundancies. In real life situations, one typically has
both: scientific computing clusters are usually composed of multiple nodes, each with
its own memory, and each node has multiple cores. It is therefore especially convenient
to use a hybrid parallelization approach, where the memory is shared within a node
and distributed across nodes.

Considering the following problem:

$$
\begin{equation*}
    \begin{cases}
    -\nabla \cdot (\mu \nabla u) + \nabla \cdot (\beta u) + \gamma u = f & \text{in } \Omega, \\
    u = g & \text{on } \Gamma_D \subset \partial \Omega, \\
    \nabla u \cdot \mathbf{n} = h & \text{on } \Gamma_N = \partial \Omega \setminus \Gamma_D
    \end{cases}
\end{equation*}
$$

Implement a solver that uses hybrid multithreading-MPI for parallelization,
relying on `deal.II`'s documentation and tutorials (with reference in
particular to tutorial steps 48 and 69).

Discuss the performance and parallel scalability of the solver, comparing to a fully
distributed MPI parallelization.

**Difficulty**: hard.
