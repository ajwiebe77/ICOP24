% whati21.m
% 
% Andrew J. Wiebe, 8 Jan 2024
% 
% Objective: Calculate total complex potential for a groundwater flow system 
%            for a wedge-shaped aquifer bounded by two constant head lines
%            based on the Holzbecher (2005) approach. Superimpose the potentials
%            related to recharge, regional flow, and well pumping.
%            Assume a permafrost (impermeable) boundary in this case.
%            Represent the Whati, NWT, peninsula and pumping well.
% 
% Notes:
%    The potential related to recharge is calculated via a numerical approach (Barba and Forsyth, 2017) that solves the Poisson equation.
%    The water table representing the regional flow field was generated using
%       cubic splines (see splinesurface3*.m)
%    Use R scripts to clean up the streamlines near the branch cut
% 
% References:
%    Barba, L.A., and Forsyth, G.F. 2017. 12 Steps to Navier-Stokes. https://nbviewer.org/github/barbagroup/CFDPython/blob/master/lessons/13_Step_10.ipynb. Cited 3 January 2024.
%    Holzbecher, E. 2005. Analytical solution for two-dimensional groundwater flow in presence of two isopotential lines. Water Resour. Res. 41, doi:10.1029/2005WR004583.
%

clear all;

tic();

%%% ---------------------------------------------------------------------------------
%%% Aquifer parameters and well locations
%%% ---------------------------------------------------------------------------------

numpts = 1000; % one less than the number of points for rz and theta; try a larger number, based on correspondence with Dr. Nagheli.

%%%% Rotation
rot_angle = -6.867509 * pi/180;
rot_matrix = [cos(rot_angle), -sin(rot_angle); sin(rot_angle), cos(rot_angle)];
x_origin = 485891; % x-location of domain origin, UTM Zone 11N
y_origin = 7001519; % y-location of domain origin, UTM Zone 11N
xw = 486009.97; % x-location of pumping well, UTM Zone 11N
yw = 7001611.41; % y-location of pumping well, UTM Zone 11N

%%%% Test parameters for Whati, NWT, using the Holzbecher (2005) approach
r_wedge = 500;
r1 = sqrt(power(y_origin - yw,2) + power(x_origin - xw, 2)); %150.42; % 100; % distance from tip approx to well point
gamma = atan((y_origin - yw) / (x_origin - xw)); %pi/6;
gamma = gamma - rot_angle; % adjust because the x-axis is rotated by rot_angle
alpha = 62.62429*pi/180; % angle between the lake boundaries, converted to radians

slope_pf = (sin(alpha) / (cos(alpha) - 1)); % slope of permafrost boundary line on 3rd side of aquifer
d = (slope_pf * r1 * cos(gamma) - slope_pf*r_wedge - r1*sin(gamma))/(sin(alpha/2) - slope_pf * cos(alpha/2));  % distance from well to pf boundary, along a line that crosses boundary at 90 degrees; only valid for alpha <= 90 deg

%%%% Image well to enforce permafrost boundary line:
%%%% Try farther away (trial and error) to account for regional gradient and obtain a stagnation point on the permafrost boundary line:
x2 = r1 * cos(gamma) + 4 * d * cos(alpha/2);
y2 = r1 * sin(gamma) + 4 * d * sin(alpha/2);

Qwell = 0.00059; % pumping rate in m3/s
b = 9; % m
K = 0.15/b; % m/s % K = T(cooper-jacobstraightline)/b

z_well = r1 * cos(gamma) + j * r1 * sin(gamma);
z_well2 = x2 + j * y2;

phi_0 = 0; % base value for total potential

%%% ---------------------------------------------------------------------------------
%%% Calculate the potential related to groundwater recharge
%%% Use numerical approach outlined by Barba and Forsyth (2017)
%%% ---------------------------------------------------------------------------------

% CASE 1: With Permafrost (PF)

nx = numpts; ny = numpts; nt = 100;
xmin = 0; xmax = r_wedge;
ymin = 0; ymax = r_wedge;

dx = (xmax-xmin)/(nx-1);
dy = (ymax-ymin)/(ny-1);

%%%% PF boundary:
% intercept_pf = xmax * (sin(alpha) - slope_pf *cos(alpha)); % intercept_pf = -slope_pf*xmax;

%initialize
p = zeros(nx,ny); % analogous to capital phi
pd = zeros(nx,ny);
x=xmin:dx:xmax;
y=ymin:dy:ymax;

%% source term (RECHARGE)
recharge_rate = 50; % 300; % mm/yr
b = - recharge_rate * 0.001 / (365 * 24 * 60 * 60); % 50 mm/yr, converted to m/s

%%% specify region for solution (i.e., not in the lake)

for nt = 1:nt
	pd = p;
	
	for i = 2:nx-1
		for ii = 2:min(floor((i * dx) * tan(alpha)/dy), ny-1)
			%%% try adding permafrost boundary:
			% if(i < (ii*dy - intercept_pf) / (slope_pf * dx))
				p(i,ii) = ((pd(i+1,ii)+pd(i-1,ii))*dy^2+(pd(i,ii-1)+pd(i,ii+1))*dx^2 - b*dx^2*dy^2)/(dx^2+dy^2)/2; % with constant b
			% end
			
			%%% Recharge extends to upland area outside unfrozen region (same as no permafrost case)
			p(i,ii) = ((pd(i+1,ii)+pd(i-1,ii))*dy^2+(pd(i,ii-1)+pd(i,ii+1))*dx^2 - b*dx^2*dy^2)/(dx^2+dy^2)/2; % with constant b
		end
	end
	
	%% change the boundary conditions - BC's for permafrost
	p(1, 2:ny-1) = p(2, 2:ny-1); % first row = second row
	p(nx, 2:ny-1) = p(nx-1, 2:ny-1); % last row = penultimate row
end

save('rechargepot_1000x1000pf.txt', 'p'); % let the recharge extend
% load('rechargepot_1000x1000pf.txt', 'p'); % let the recharge extend past the permafrost boundary to avoid edge effects

%%% ---------------------------------------------------------------------------------
%%% Construct regional gradient from water levels developed in splinesurface3*.m
%%% ---------------------------------------------------------------------------------

phi_z = zeros(numpts,numpts);

load('reg_background_wls_1000zzz_00001.txt', 'zzz'); % now smoothed using subsamplint and interp2() function

% Note: use equation developed from Strack (2017), and discretize: capPhi = 0.5*K*h*h; Qx0 = -dai(capPhi)/dai(x) = -d(capPhi)/dx = (capPhi2 - capPhi1) / dx
regflowx = zeros(numpts,numpts);
regflowy = zeros(numpts,numpts);

for i = 2:numpts
	for ii = 2:numpts
		regflowx(ii,i) = - K * 0.5 * (power(zzz(ii,i),2) - power(zzz(ii,i-1),2)) / (r_wedge / numpts);
		regflowy(ii,i) = - K * 0.5 * (power(zzz(ii,i),2) - power(zzz(ii-1,i),2)) / (r_wedge / numpts);
	end
end

regflowx(:,1) = regflowx(:,2);
regflowy(:,1) = regflowy(:,2);
regflowx(1,:) = regflowx(2,:);
regflowy(1,:) = regflowy(2,:);

%%% ---------------------------------------------------------------------------------
%%% Calculate the total potential by superposition terms representing the recharge, regional gradient, and wells
%%% Equation adapted from Holzbecher (2005) Equation 8, noting from Nagheli et al. (2020) that the exponent on the complex coordinates is the aquifer angle at the origin divided by pi
%%% ---------------------------------------------------------------------------------

for i = 1:numpts
	for ii = 1:numpts
		z_x(ii,i) = r_wedge * (i - 1) / (numpts - 1);
		z_y(ii,i) = r_wedge * (ii - 1) / (numpts - 1);
		z = z_x(ii,i) + z_y(ii,i) * j;
		
		% NOTE: restrict regional gradient to the aquifer region (i.e., do not apply within the lake)
		
		if(z_y(ii,i) <= z_x(ii,i) * tan(alpha)) % within aquifer
			phi_z(ii,i) = phi_0 + p(ii,i) - (regflowx(ii,i) - j * regflowy(ii,i))*z + (Qwell / (2 * pi)) * (log(power(z,pi/alpha) - power(z_well,pi/alpha)) - log(power(z,pi/alpha) - power(conj(z_well),pi/alpha)) + log(power(z,pi/alpha) - power(z_well2,pi/alpha)) - log(power(z,pi/alpha) - power(conj(z_well2),pi/alpha)));
			
		else % within the lake - remove regional gradient term and recharge term
			phi_z(ii,i) = phi_0 + (Qwell / (2 * pi)) * (log(power(z,pi/alpha) - power(z_well,pi/alpha)) - log(power(z,pi/alpha) - power(conj(z_well),pi/alpha)) + log(power(z,pi/alpha) - power(z_well2,pi/alpha)) - log(power(z,pi/alpha) - power(conj(z_well2),pi/alpha)));
		end
		
	end
end

phi_z(1,:) = phi_z(2,:); % This adjusts for anomalous first row

phi_z_R = real(phi_z);
phi_z_I = imag(phi_z);
save('phi_z1000x1000real_pf.csv', 'phi_z_R');
save('phi_z1000x1000imag_pf.csv', 'phi_z_I');

% figure; surf(z_x,z_y,p);
% figure; surf(z_x, z_y, phi_z_R);

%%%% FIND STAGNATION POINT --------------------------------------------------------------
%%%% With help from the GNU Octave manual page at: https://octave.org/doc/v5.2.0/Solvers.html

maxrealz = 30; % realizations
maxIter = 250; % iterations per realization

xf = zeros(maxIter,1);

gradx = zeros(numpts,numpts);
grady = zeros(numpts,numpts);

for i = 1:(numpts-1)
	for ii = 1:(numpts-1)
		gradx(i,ii) = real(phi_z(i,ii+1)) - real(phi_z(i,ii));
		grady(i,ii) = real(phi_z(i+1,ii)) - real(phi_z(i,ii));
	end
end

gradx(:,numpts) = gradx(:,numpts-1);
grady(:,numpts) = grady(:,numpts-1);
gradx(numpts,:) = gradx(numpts-1,:);
grady(numpts,:) = grady(numpts-1,:);


%%% could possibly use the gradient function here: [dx, dy] = gradient(real(phi_z))
%%% not yet tested: % [gradx, grady] = gradient(real(phi_z));

% https://stackoverflow.com/questions/13136987/octave-compute-gradient-of-a-multi-dimensional-function-at-a-particaular-point
gradient3 = @(i,ii) [gradx(i,ii), grady(i,ii)]; % a wrapper function to allow both the x and y direction gradients to be returned for a pair of x and y coordinates

lowestgrad = [0,0,100]; % [x,y,dummy value]

for realz = 1:maxrealz

	x = zeros(maxIter,2);

	xy = [rand() * r_wedge, rand() * r_wedge]; % exact (x,y) location
	x(1,1:2) = [ceil(numpts * xy(1,2) / r_wedge), ceil(numpts * xy(1,1) / r_wedge)]; % indices of cell for exact location
	
	for i = 2:maxIter
		% https://www.khanacademy.org/math/multivariable-calculus/applications-of-multivariable-derivatives/optimizing-multivariable-functions/a/what-is-gradient-descent
		
		g = gradient3(x(i-1,1),x(i-1,2));
		mult = min(abs(5E-1 / g(1,1)), abs(5E-1 / g(1,2)));
		
		xy = xy(1,1:2) .+ mult * g;		% update exact location
		row = min(1000, ceil(numpts * xy(1,2) / r_wedge));
		row = max(1, row);
		col = min(1000, ceil(numpts * xy(1,1) / r_wedge));
		col = max(1,col);
		
		x(i,1:2) = [col, row];

	end

	dphi = gradient3(x(maxIter,1), x(maxIter,2));
	resultant3 = sqrt(power(dphi(1,1),2) + power(dphi(1,2),2));

	% ignore results very close to the apex of the aquifer, where the water table gradient is very slow by design
	if (resultant3 < lowestgrad(1,3)) && (x(maxIter,1) > 10) && (x(maxIter,2) > 10) % choose a point that's not near the peninsula tip
		lowestgrad(1,1) = x(maxIter,1);
		lowestgrad(1,2) = x(maxIter,2);
		lowestgrad(1,3) = resultant3;
	end
	
end

lowestgrad

xy_stg = [z_x(lowestgrad(1,1), lowestgrad(1,2)), z_y(lowestgrad(1,1), lowestgrad(1,2))]

save xy_stg_pf.txt xy_stg;

elapsed_time = toc()