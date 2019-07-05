% Sat  8 Apr 20:57:47 CEST 2017
%
%% analytic solution to the river tide formulated as boundary value problem
%% in a river with finite length
%%
%% c.f. Godin 1986
%%
% TODO make velocity scale not constant
% TODO make friction non-linear
% TODO solve for two species that interact
% TODO allow for width and depth variation
function [x, z] = damped_wave_bvp(L)
	
	g = 9.81;
	T = 86400;
	o = 2*pi/T;
	cD = 5e-4;
	H = 10;
	z0 = 1;
	u0 = 1;
	% damping coefficient according to lorentz
	r = -cD*(8/(3*pi))*u0;
	% stabilistation
	s = 0;
	
	a = 1/(H*g)*(-o^2 + r/H*1i*o)	
	
	% boundary value problem
	% (not that godin has z(L) = z and z(0) = 0, leading to a mirroredt expression of the solution
	% z(0) = z0 -> 
	ac     = z0;
	% z(L) = 0
	sa     = sqrt(-a);
	as     = -ac*cos(sa*L)/sin(sa*L);
	z      = as*sin(sa*x) + ac*cos(sa*x);
	z      = z0*(-sin(sa*x)*cos(sa*L)/sin(sa*L) + cos(sa*x));
	z(:,2) = z0*(sin(sa*(L-x))/sin(sa*L));
	%z(:,2) = z0*(sin(sa*(-x))/sin(sa*L)); % L>>x
	dz_dx = (z(2,1)-z(1,1))/(x(2)-x(1))
	
	% Godin z(L) = z0*1, z(0) = 0
	ac = 0;
	as = 1/sin(sqrt(a)*L);
	z  = z0*(as*sin(sa*x) + ac*cos(sa*x));
	z  = z0*(sin(sa*x)/sin(sa*L));
	
	

	% solution to BVP	
	a  = a;
	sa = sqrt(-a);
	
	x  = linspace(0,L,n)';
	% cos part
	ac = z0;
	% sin part
	as = -ac*cos(sa*L)/sin(sa*L);
	z  =  as*sin(sa*x) + ac*cos(sa*x);
	
end % damped_wave_bvp


