README.txt

Author: Andrew J. Wiebe
Date: 10 Jan 2024

ICOP24 Repository (ajwiebe77/ICOP24)

Code for GNU Octave 6.4.0 (Eaton et al., 2020)
(Note: This code (.m files) is likely nearly compatible with MATLAB; it may require minor updates. This has not been tested.)
and code for R (R Core Team, 2022)


Introduction:

The code in this repository (ajwiebe77/ICOP24) contains calculations related to analytical/numerical simulations of groundwater capture zones for a well on a wedge-shaped aquifer bounded by two or three equipotential lines (two lake boundaries and either a permafrost boundary or a semi-infinite boundary). The two cases addressed are: 1) A permafrost boundary at x = 500 m along the x-axis, and 2) No permafrost boundary.

The idea for this work is based on a paper by Nagheli et al. (2020), but the method was coded with alternative (dimensioned) equations in order to more easily include recharge into the superposition calculation. It might be possible to derive an analytical expression for the recharge for use with the Nagheli et al. (2020) equations, but this was not attempted.

The regional gradient is first constructed (splinesurface*.m) from a water table elevation matrix constructed using cubic splines. A regional background water table gradient is assumed (at a radius of 1000 m from the origin) and a gradient of 0 m/m is assumed for flow lines perpendicular to and crossing the lake boundaries. The text file generated is then read by one of the two .m scripts (one for each case noted above). The script (whati21.m for 'permafrost' or whati20.m for 'no permafrost') calculates or reads a potential field (matrix of values) related to groundwater recharge, reads the regional water table matrix file and calculates the x- and y-direction gradients, then superimposes the potentials related to recharge, regional flow, and well pumping. Finally, a gradient ascent method is used to approximate the location of a stagnation point (only zero or one is assumed for this system).

The field site for the application of this method is the community of Whatì, NWT, Canada (117’16”19.9 °W, 63’8”38.3 °N). The well pumping rate and aquifer thickness were estimated based on Stanley Associates (1987). The distance to permafrost was loosely based on geotech studies by Thurber (1981) and Hardy (1991).

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
- Ensuring visual continuity of streamlines across a branch cut was performed via the sequential contouring method by Holzbecher (2018). This method was previously coded (Wiebe and McKenzie (2022)) - see the ajwiebe77/YWVC respository - and was modified slightly for use here.


Known Limitations:

- The background regional water table has some issues outside of a radius of 1000 m near the edge of the 1000 m by 1000 m area. This may be adjusted by ghosting in some points along the outer edges. This would also be helpful along the two boundary arms of the aquifer wedge for distances from the origin of > 500 m. A more realistic water table would slope down those lines and then reach an elevation of zero for the 0 to 500 m distance along each boundary line.
- Future work could attempt to more reasonably follow the lake shoreline and alter the geometry away from a simple wedge. The Nagheli et al. (2020) approach for irregular polygon-shaped aquifers may be useful in association with this (recharge and regional gradient terms would need to be added to their analytical solutions). It might be possible to identify stagnation points and therefore capture zone boundaries if the wider flow field was considered. (Regional flow may branch with some flowing underneath the community and toward the well and some flowing directly into the lake on either side of the peninsula without entering the peninsula wedge).
- The approach of using an image well fails to create a no-flow/impermeable permafrost boundary in Case (1). This is likely due to the presence of a regional flow field, because Mahdavi (2021), for instance, illustrates that an impermeable boundary can be created for a wedge aquifer using an image well in the absence of a regional gradient.
- The gradient ascent method (inspired by a gradient descent method; e.g., Khan Academy [2024]) is far from perfect and could be improved. Problems: 1) the method only identifies one point, which may be erroneous and located on or very close to a boundary or edge of the domain; 2) there is no feedback to suggest whether the point is likely erroneous or a good estimate; 3) the step size may not be optimal; 4) the algorithm could pass through a stagnation point and keep going rather than identifying where best to stop; and 5) the approach of using a random initial location may not be optimal.


References:

Barba, L.A., and Forsyth, G.F. 2017. 12 Steps to Navier-Stokes. https://nbviewer.org/github/barbagroup/CFDPython/blob/master/lessons/13_Step_10.ipynb. Cited 3 January 2024.

Eaton, J.W., Bateman, D., Hauberg, S., and Wehbring, R., 2020. GNU Octave version 6.4.0 manual: A high-level interactive language for numerical computations. https://docs.octave.org/v6.4.0/index.html. Cited 3 Jan 2024.

Hardy BBT Ltd. (Hardy). 1991. Geotechnical investigation: Proposed Community Office, Lac La Martre, NWT. Prepared for: GNWT Dept. of Public Works, Architectural Division, Yellowknife, NWT. Hardy BBT Ltd.: Yellowknife, NWT. Mar 1991.

Holzbecher, E. 2005. Analytical solution for two-dimensional groundwater flow in presence of two isopotential lines. Water Resour. Res. 41, doi:10.1029/2005WR004583.

Holzbecher, E. 2018. Streamline visualization of potential flow with branch cuts, with applications to groundwater. J. Flow Vis. Image Process. 25(2), 119–144. doi:10.1615/JFlowVisImageProc.2018025918.

Khan Academy, 2024. Gradient descent. https://www.khanacademy.org/math/multivariable-calculus/applications-of-multivariable-derivatives/optimizing-multivariable-functions/a/what-is-gradient-descent. Cited 10 Jan 2024.

Mahdavi, A. 2021. Response of dual-zone heterogeneous wedge-shaped aquifers under steady-state pumping and regional flow. Advan. Water Resour. 147, 103823, doi:10.1016/j.advwatres.2020.103823.

Nagheli, S., Samani, N.,and Barry, D.A., 2020. Capture zone models of a multi-well system in aquifers bounded with regular and irregular inflow boundaries. J. Hydrol. X 7, 100053. doi:10.1016/j.hydroa.2020.100053.

R Core Team. 2022. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/. Cited 17 March 2022.

Stanley Associates Engineering Ltd. (Stanley Associates). 1987. Lac La Martre - Well Construction, Draft Report. Prepared for Government of Northwest Territories, Engineering Division, Department of Public Works and Highways. Stanley Associates Engineering, Ltd.: Yellowknife, NT, Canada. 23 Mar 1987.

Strack O.D.L. 1989. Groundwater Mechanics. Prentice-Hall, Englewood Cliffs, NJ, USA.

Strack, O.D.L. 2017. Analytical groundwater mechanics. Cambridge Univerisity Press, Cambridge, United Kingdom. doi:10.1017/9781316563144.

Thurber Consultants Ltd. (Thurber). 1981. Lac La Martre School: Geotechnical Investigation. Prepared for: Dept. of Public Works, Government of the Northwest Territories. Thurber Consultants Ltd.: Edmonton, AB. 9 March 1981.

Van Loan, C.F. 2000. Introduction to Scientific Computing. 2nd Ed, Prentice-Hall, Upper Saddle River, NJ, USA.

Wiebe, A.J., and McKenzie, J.M. 2022. An Open-Source Web Tool for Visualizing Estimates of Well Capture Zones Near Surface Water Features. Poster presentation at: AGU Fall Meeting 2022, Chicago, IL, USA, 12-16 Dec 2022. doi:10.22541/essoar.167267811.10671930/v1.

