% Tue 19 Jul 17:27:31 CEST 2016
%% phase lag of high and low water
%%
%% phi : u_river/u_tide < 1
%%
%% delta_eps_hw = omega*(t_hws - t_hw)
%% delta_eps_hw = omega*(t_lws - t_lw)
%%
%% c.f. savenije
function [delta_eps_hw, delta_eps_lw] = savenije_phase_lag(T,h,U_t,U_r,a,delta)
	g = 9.81;

	phi = U_r./U_t;
	% v : tidal velocity amplitude
	
	% angular frequency of main constituent
	omega = 2*pi./T;

	% classical wave celertiy
	c0 = sqrt(g*h);

	% TODO
	c = c0;
	
	% phi 		phase lag angle with respect to the reference station
	% lambda	celerity number (3.22,3.52)
	lambda = c0./c;

	% gamma		shape number (3.53)
	gamma = c0./(omega*a);

	% delta		damping number
	% F		Froude number	(1.3)
	% does river discharge needs to be considered here? Savenije leaves this unclear
	F = U_t./c0;
	
	% without 
%	if (phi > 1)
%		delta_eps_lw = pi/2;
%		delta_eps_hw = -pi/2;
%	else
		delta_eps_hw = -asin(phi) + atan(lambda./(gamma - delta).*(1 + gamma./lambda.*F.*phi*pi - 2*delta./lambda));
		delta_eps_lw = +asin(phi) + atan(lambda./(gamma - delta).*(1 + 2*gamma./lambda.*F.*phi*pi));
%	end
end	





