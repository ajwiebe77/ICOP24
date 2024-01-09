README.txt

Author: Andrew J. Wiebe
Date: 3 Jan 2024

ICOP24 Repository

Code for GNU Octave 6.4.0 (Eaton et al., 2020)
(Note: The code is likely nearly compatible with MATLAB; may require minor updates. This has not been tested.)
and code for R (R Core Team, 2022)


Introduction:

The code in this repository (ajwiebe77/ICOP24) contains calculations related to analytical/numerical simulations of groundwater capture zones for a well on a wedge-shaped aquifer bounded by two or three equipotential lines (two lake boundaries and either a permafrost boundary or a semi-infinite boundary). The two cases addressed are: 1) A permafrost boundary at x = 500 m along the x-axis, and 2) No permafrost boundary.

The idea for this work is based on a paper by Nagheli et al. (2020), but the method was coded with alternative (dimensioned) equations in order to more easily include recharge into the superposition calculation. It might be possible to derive an analytical expression for the recharge for use with the Nagheli et al. (2020) equations, but this was not attempted.

The regional gradient is first constructed (splinesurface*.m) from a water table elevation matrix constructed using cubic splines. A regional background water table gradient is assumed (at a radius of 1000 m from the origin) and a gradient of 0 m/m is assumed for flow lines perpendicular to and crossing the lake boundaries. The text file generated is then read by one of the two .m scripts (one for each case noted above). The script (whati20.m for no permafrost or whati21.m for permafrost) calculates or reads a potential field (matrix of values) related to groundwater recharge, reads the regional water table matrix file and calculates the x- and y-direction gradients 



Highlights:

Users may find the following aspects of this project useful:
- The GNU Octave code applies the superposition method to several potential fields (related to groundwater recharge, regional unconfined water table gradient, wells) and generates the resulting complex potential field.
- The R scripts whati9pf.R and whati9nopf.R contain code that saves contour lines data to .shp files (as points). The 'order' and 'group' fields allow the points to be processed in GIS software (e.g., 'Points to Path' operation in QGIS 3.16.9) to create polylines. Some editing of the resulting polylines may be necessary.


ICOP24 Coding Notes:

- Regional gradient (beta) angles in the GNU Octave code are in radians counterclockwise from east, rather than as the values stated clockwise from north in the paper.
- The header lines and first (blank) column need to be removed from the csv files prior to reading them using the appropriate R script (whati9pf.R or whati9nopf.R)
- Regional gradient potential term in the superposition equation - from Strack (2017) and Strack (1989)
- Numerical solution method for solving Poisson equation for recharge is from Barba and Forsyth (2017)
- The superposition equation (for complex potential) for an aquifer with two isopotential/equipotential lines was adapted from Holzbecher (2005), with thanks to Nagheli et al. (2020) for the understanding of the powers of the complex numbers (pi/alpha) in the logarithm functions.
- Ensuring visual continuity of streamlines across a branch cut was performed via the sequential contouring method by Holzbecher (2018). This method was previously coded by Wiebe and McKenzie (2022) - see the ajwiebe77/YWVC respository - and was 
modified slightly for use here.



References:

Barba, L.A., and Forsyth, G.F. 2017. 12 Steps to Navier-Stokes. https://nbviewer.org/github/barbagroup/CFDPython/blob/master/lessons/13_Step_10.ipynb. Cited 3 January 2024.

Eaton, J.W., Bateman, D., Hauberg, S., and Wehbring, R., 2020. GNU Octave version 6.4.0 manual: A high-level interactive language for numerical computations. https://docs.octave.org/v6.4.0/index.html. Cited 3 Jan 2024.

Holzbecher, E. 2005. Analytical solution for two-dimensional groundwater flow in presence of two isopotential lines. Water Resour. Res. 41, doi:10.1029/2005WR004583.

Holzbecher, E. 2018. Streamline visualization of potential flow with branch cuts, with applications to groundwater. J. Flow Vis. Image Process. 25(2), 119–144. doi:10.1615/JFlowVisImageProc.2018025918.

Nagheli, S., Samani, N., and Barry, D.A., 2020. Capture zone models of a multi-well system in aquifers bounded with regular and irregular inflow boundaries. J. Hydrol. X 7, 100053. doi:10.1016/j.hydroa.2020.100053.

R Core Team. 2022. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/. Cited 17 March 2022.

Strack O.D.L. 1989. Groundwater Mechanics. Prentice-Hall, Englewood Cliffs, NJ, USA.

Strack, O.D.L. 2017. Analytical groundwater mechanics. Cambridge Univerisity Press, Cambridge, United Kingdom. doi:10.1017/9781316563144.

Van Loan, C.F. 2000. Introduction to Scientific Computing. 2nd Ed, Prentice-Hall, Upper Saddle River, NJ, USA.

Wiebe, A.J., and McKenzie, J.M. 2022. An Open-Source Web Tool for Visualizing Estimates of Well Capture Zones Near Surface Water Features. Poster presentation at: AGU Fall Meeting 2022, Chicago, IL, USA, 12-16 Dec 2022. doi:10.22541/essoar.167267811.10671930/v1.

