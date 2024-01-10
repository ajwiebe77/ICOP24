% splinesurface3_dhdl00001.m
% 
% Andrew J. Wiebe, 6 Jan 2024
% 
% Objective: Use cubic splines to generate a surface representing the water table within the Whati, NWT, peninsula
% 
% Notes:
%    Use spline math equations from Van Loan (2000)
%    Subsample a grid of points from the results of the griddata function (linear) and use interp2 with cubic interpolation to get a smoother version of the water table.
%    There may be some problems in the contour map along the two lake boundaries
%    zz_list = (x, y, water table elevation)
%    
% References:
%    Van Loan, C.F. 2000. Introduction to Scientific Computing. 2nd Ed, Prentice-Hall, Upper Saddle River, NJ, USA.
% 

clear all;

tic();

%%% No PF case

% -----------------------------------------------------------------------

%%% Equations of curves starting at various fractions between centre-point along the arc and each of the lake boundary endpoints

%%% Upper curve
r = 1000; % radius (larger than than the desired calculation area, which has a radius of 500 m)
n = 2000; % number of points % want 1000 points between x = 0 and 500 m
alpha  = 62.62429*pi/180;
beta = pi + alpha/2;

outerarc = [0.5 * r * cos(0:alpha/10:alpha);  0.5 * r * sin(0:alpha/10:alpha)]';

% figure;
% hold on;
% plot([r*cos(alpha), 0, r], [r*sin(alpha), 0, 0], 'k-', outerarc(:,1), outerarc(:,2), 'k--', [0, r*cos(alpha/2)], [0, r*sin(alpha/2)]);
% axis([0,r,0,r]);

counter = 1;

% Centreline -------------------------------------------------
y1c = 0 ; % S1(x1)
dh_dl = 1E-5;
x1c = 0;
x2c = r; %r*sqrt(2); % 500;
deltaX1c = x2c - x1c; % goes from tip to midpoint of arc

%%% assume the overall gradient between the peninsula tip and the midpoint along the outer arc is equal to the regional gradient
y1_primec = dh_dl;

y2c = y1_primec * deltaX1c + y1c;

g1c = 0;
g2c = dh_dl;
s1c = g1c;
s2c = g2c;

a1c = y1c;
b1c = s1c;
c1c = (y1_primec - s1c) / deltaX1c;
d1c = (s2c + s1c - 2 * y1_primec) / power(deltaX1c,2);

stepc = (1/n) * r; % 1;
xc = 0:stepc:deltaX1c;
maxzz = 0;

for i = 1:length(xc)
	S1c(i) = a1c + b1c * (xc(i)) + c1c * power(xc(i),2) + d1c*power(xc(i),2) * (xc(i) - x2c);
	
	zz_list(counter,1:3) = [xc(i) * cos(alpha/2), xc(i) * sin(alpha/2), S1c(i)];
	counter = counter + 1;
	
	if S1c(i) > maxzz
		maxzz = S1c(i);
	end
end

% plot(zz_list(:,1),zz_list(:,2),'.');
% plot3(zz_list(:,1),zz_list(:,2),zz_list(:,3));

% Lake Boundaries -------------------------------------------------

xstep = n/10;

for ii = 1:xstep:r
	zz_list(counter,1:3) = [ii, 0, 0];
	counter = counter + 1;
end

for ii = 1:xstep:(r*cos(alpha))
	zz_list(counter,1:3) = [ii, ii * tan(alpha), 0];
	counter = counter + 1;
end


%%% Splines along outer extent (infinite side) -------------------------------------------------
%%% Upper (mirror image across centreline for lower)
y1o = 0 ; % S1(x1)
x1o = 0;
x2o = r*alpha / 2;
deltaX1o = x2o - x1o; % goes from tip to midpoint of arc

%%% Assume the overall gradient between the peninsula tip and the midpoint along the outer arc is equal to the regional gradient
y1_primeo = (maxzz - y1o) / deltaX1o;
y2o = maxzz;

s1o = 0;
s2o = 0;

a1o = y1o;
b1o = s1o;
c1o = (y1_primeo - s1o) / deltaX1o;
d1o = (s2o + s1o - 2 * y1_primeo) / power(deltaX1o,2);

stepo = deltaX1o / ((((31/32) - (22/32)) / (1/32)) + 1);
xo = stepo:stepo:deltaX1o;

for i = 1:length(xo)
	S1o(i) = a1o + b1o * (xo(i)) + c1o * power(xo(i),2) + d1o*power(xo(i),2) * (xo(i) - x2o);
end

% figure; plot(xo,S1o);



% Arcs --------------------------------------------------------
% Use a spline to determine the (x,y) coordinates of points along each path line

kstep = 1/32;

for k = (22/32):kstep:(31/32)
	x2 = r * cos(k*alpha);
	y2 = r * sin(k*alpha);

	x1 = (k - 0.5)*r*cos(alpha);
	y1 = (k - 0.5)*r*sin(alpha);

	slopebndy2 = tan(alpha);

	s1 = - 1/tan(alpha);
	s2 = tan(beta);

	deltaX1 = x2 - x1;

	y1prime = (y2 - y1) / deltaX1;

	a1 = y1;
	b1 = s1;
	c1 = (y1prime - s1)/deltaX1;
	d1 = (s2 + s1 - 2 * y1prime) / power(deltaX1,2);
	
	step = (1/n) * r; 1; %0.1;
	x = zeros(1, floor((x2 - x1) / step) + 1);
	S1 = zeros(1, length(x));
	x_lower = zeros(1, length(x));
	y_lower = zeros(1, length(x));
	x = x1:step:x2;

	for i = 1:length(x)

		S1(i) = a1 + b1 * (x(i) - x1) + c1 * power(x(i) - x1,2) + d1*power(x(i) - x1,2) * (x(i) - x2);
		
	end
	
	x_upper = x;
	y_upper = S1;
	
	%%%% Lower curve

	slope2 = - 1 / (tan(alpha/2));

	for i = 1:length(x)
		bi = S1(i) - slope2 * x(i);
		x_cl = -bi / (slope2 - tan(alpha / 2));
		y_cl = x_cl * tan(alpha/2);
		d = sqrt(power(y_cl - S1(i),2) + power(x_cl - x(i),2));
		aa = power(slope2,2) + 1;
		bb = 2 * slope2 * bi - 2 * slope2 * S1(i) - 2 * x(i);
		cc = power(bi,2) + power(S1(i),2) - 2 * bi *S1(i) + power(x(i),2) - 4 * power(d,2);
		x_lower(i) = (-bb + sqrt(power(bb,2) - 4 *aa * cc)) / (2 * aa);
		y_lower(i) = slope2 * x_lower(i) + bi;
	end 
	
	% outerarc = [r * cos(0:alpha/10:alpha);  r*sin(0:alpha/10:alpha)]';
	
	% plot(x_upper,y_upper);
	% plot(x_lower,y_lower);
	
	% -------------------------------------------------------------
	%%% Use a spline to estimate the z values at the points along each path line
	
	y1h = 0; % S1(x1)
	x1h = 0;
	x2h = 0;
	
	for i = 2:length(x_upper)
		x2h = x2h + sqrt(power(x_upper(i) - x_upper(i-1),2) + power(y_upper(i) - y_upper(i-1),2));
	end
	
	deltaX1h = x2h - x1h;
	
	%%% Assume the overall gradient between the peninsula tip and the midpoint along the outer arc is equal to the regional gradient
	y1 = 0;
	y2 = S1o(floor(length(S1o) - (k - 22/32) / kstep)); % use outer spline; upper
	y1_primeh = (y2 - y1) / deltaX1h;
	
	s1h = 0;
	s2h = y1_primeh;

	a1h = y1h;
	b1h = s1h;
	c1h = (y1_primeh - s1h) / deltaX1h;
	d1h = (s2h + s1h - 2 * y1_primeh) / power(deltaX1h,2);
	
	steph = deltaX1h / (length(x_lower) - 1); %1;
	
	xh = zeros(1, floor(deltaX1h / steph)  + 1);
	S1h = zeros(1, length(xh));
	
	xh = 0:steph:deltaX1h;
	
	for i = 1:length(xh)
		
		S1h(i) = a1h + b1h * (xh(i)) + c1h * power(xh(i),2) + d1h*power(xh(i),2) * (xh(i) - x2h);
		
		zz_list(counter,1:3) = [x_lower(i), y_lower(i), S1h(i)];
		counter = counter + 1;
		
		zz_list(counter,1:3) = [x_upper(i), y_upper(i), S1h(i)];
		counter = counter + 1;
		
	end
	
	% figure; plot(xh,S1h);
	% axis([0, n, 0, 0.01]);
	% title(num2str(k));
	
	
end


%%% Interpolate via griddata (a linear function; cubic interpretation not currently supported):

xi = linspace(0, r/2, n/2);
yi = linspace(0, r/2, n/2);
[xxx, yyy] = meshgrid(xi,yi);
zzz = griddata(zz_list(:,1), zz_list(:,2), zz_list(:,3), xxx, yyy); % problem: griddata does not have the cubic interpolation implemented

%%% Subsample the results and then use interp2 to generate a smoother surface

zzz(isnan(zzz)) = 0;

rowsub = [1,200:200:1000];
colsub = rowsub;

xsub = xxx(1, colsub);
ysub = yyy(rowsub,1);
zsub = zzz(rowsub, colsub);

xi2 = xi(1,:);
yi2 = yi(1,:)'; % transpose
zzz2 = interp2(xsub,ysub,zsub,xi2,yi2,'cubic');

[xmat, ymat] = meshgrid(xsub,ysub);
% mesh(xi2,yi2,zzz2);
% hold on; plot3(xmat,ymat,zsub,'o');

figure; contour(xi2,yi2,zzz2);

zzz = zzz2;

save('reg_background_wls_1000zzz_00001.txt', 'zzz'); % interpolated values

% surf(xxx,yyy,zzz)

elapsed_time = toc()
