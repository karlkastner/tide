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
function [x, z,dz_dx] = river_tide_godin(L,z0,u0,H,cD,T)
	
	g = Constant.gravity;
	uif (nargin()< 1 || isempty(L))
		L = 5e5;
	end
	if (nargin()< 2 || isempty(z0))
		z0 = 1;
	end
	if (nargi()< 3 || isempty(u0))
		u0 = 1;
	end
	if (nargin()< 4 || isempty(H))
		% default depth : 10m
		H = 10;
	end
	if (nargin < 5 || isempty(cD))
		cD = 2.5e-2; %5e-4;
	end
	if (nargin()<6 || isempty(T))
		% default: semidiurnal tida
		T = 0.5*Constant.SECONDS_PER_DAY;
	end
	n = 100;

	o = 2*pi/T;

	% damping coefficient according to lorentz
	r = -cD*(8/(3*pi))*u0;

	% stabilistation
	s = 0;
	
	a = 1/(H*g)*(-o^2 + r/H*1i*o)	
	% (not that godin has z(L) = z and z(0) = 0, leading to a mirrored expression of the solution
	
	% solve boundary value problem
	% z(0) = z0, z(L) = 0
	[x,z,dz_dx] = bvp_1d(z0,L,a)
	% z = bvp_1d(z0,L,a)
end % damped_wave_bvp

