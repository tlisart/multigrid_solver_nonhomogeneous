This is a project made for the course MATH-H401 in the first year master degree
at the Ecole Polytechnique de Bruxelles, the main goal of this project is to
build a multi-grid solver for linear symmetric sparse problems, resulting from
a discretization.

Before compiling make sure you have the following libraries placed at the root
folder containing "main.c" :

  -PRIMME  (pre-compiled)
  -UMFPACK (pre-compiled)

To compile : "make" command, the executable is by default "./main"
To clean all folders : "make clean" command

PRIMME will not run as long as the python script displaying the residual evolution is not closed

(Projet 9)
================================================================================
Content :
  - General multi-grid solver adapted to the project
  - PRIMME resolution with preconditionner acceleration
================================================================================

Use :
  - All modifiable variables are at the top of the main.c file
    - Changing to m_depth = 9 for multigrid solver (called in PRIMME)

  - General parameters :
    - ptype             --> changing smoothing method
    - relaxationOn      --> Calculate the best \tau for optimal relaxation
    - stabCheck         --> checks the stability criterion after computing
    - primmePrecond     --> Skips or not the eigenvalue problem
    - primmePrecondType --> Choose what preconditionner to use in PRIMME

  - Graphical options :
    - showDiffusionMap       --> Heat map diffusion problem after computing
    - showDiffusionSurface   --> Shows the surface of the diffusion problem
    - showEigenvalues        --> Displays the first mode given by the PRIMME eigenvalue solver
    - showResGraph           --> Makes a graph of the residual norm
================================================================================
