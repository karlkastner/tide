% Sat  8 Apr 20:57:47 CEST 2017
%% solved damped wave equation
%% z'' + a z = 0
%% z(0) = z0, z(L) = 0
function [x,z,dz_dx] = bvp_1d(z0,L,a,n)
		x      = linspace(0,L,n)';
		ac     = z0;
		sa     = sqrt(-a);
		as     = -ac*cos(sa*L)/sin(sa*L);
		z      =  as*sin(sa*x) + ac*cos(sa*x);
		z(:,2)  =  z0*(-sin(sa*x)*cos(sa*L)/sin(sa*L) + cos(sa*x));
		z(:,3) = z0*(sin(sa*(L-x))/sin(sa*L));
		%z(:,2) = z0*(sin(sa*(-x))/sin(sa*L)); % L>>x
		dz_dx = (z(2,1)-z(1,1))/(x(2)-x(1))
end
% swap left and right
	% Godin z(L) = z0*1, z(0) = 0
	function z = bvp_1d_(z0,L,a,n)
		x  = linspace(0,L,n)';
		ac = 0;
		as = 1/sin(sqrt(a)*L);
		z  = z0*(as*sin(sa*x) + ac*cos(sa*x));
		z  = z0*(sin(sa*x)/sin(sa*L));
	
		% solution to BVP	
		%a  = a;
		%sa = sqrt(-a);
		
		% cos part
		ac = z0;
		% sin part
		as = -ac*cos(sa*L)/sin(sa*L);
		% sum
		z  =  as*sin(sa*x) + ac*cos(sa*x);
	end

