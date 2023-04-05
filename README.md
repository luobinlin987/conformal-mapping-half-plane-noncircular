# Codes in this repo are for paper in arXiv DOI:10.48550/arXiv.2209.05752

#This is the README file for the FORTRAN codes related to the manuscript "On solution of conformal mapping for a lower half plane containing a symmetrical noncircular cavity".

# The FORTRAN codes are full realization of the pseudocodes in the manuscript.

# The FORTRAN codes for Section 3 are:

1- non-convex-verification.f90

2- parameters.f90

# The FORTRAN codes for Section 4 are separated into 8 .f90 files：

1- main.f90 

2- parameters.f90 (shared with Section 3)

3- PS-1.f90 

4- PS-2.f90

5- PS-3.f90

6- PS-4.f90

7- PS-final.f90

8- Subroutine-core.f90

# In these .f90 files, files 1, 2, 3, and 8 are detailedly marked with explanatory notes to facilitate readers and users. The explanatory notes of the variables in the manuscript are written in the form of LaTeX for identification.

# The file "iteration.gif" is the data visualization of rectangular coordinate convergence, as the iterations proceed (0-3000 iterations).

# Since the solution is a modification of Particle Swarm Method, randomness is a part of the solution. Thus, there will be slightly different when running the FORTRAN program or realizing the pseudocodes with another language.
