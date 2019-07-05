% Fri 17 Nov 09:26:23 CET 2017
% Karl Kastner, Berlin
%
%% damping modulus of the tide without river flow
%% c.f. friedrichs, ippen harleman
%
% function [r r_low r_high] = damping_modulus_tide(omega,cd,h0,az1)
%% output :
%% k : wave number
%% re(k) : rate of phase change
%% im(k) : damping rate
function [k, r_low, r_high] = damping_modulus_tide(omega,cd,h0,az1)
	full = true;
	if (~issym(az1))
		g  = 9.81; % is eliminated
		Pi = pi;
	else
		syms g Pi
	end
	% celerity of the undamped wave in shallow water
%	c  = sqrt(g*h0);
	% wave number of the undamped wave
%	k  = omega/c;
	% velocity of the undamped wave
%	u1 = az1*sqrt(g/h0);
	% discharge of the undamped wave
%	qt = u1*h0;
	qt = az1*sqrt(g*h0);
 
	% This is from friedrichs, but might not be correct:
%	r_ = 8/(3*Pi)*cd*u1/h0;
%	r  = -(k.*(omega./r_-sqrt(omega.^2./r_.^2+1)));

	re = -omega.^2./(g.*h0);
	im = 8/(3*Pi).*omega.*cd./g.*qt./h0.^3;
	if (~issym(re))
		k = sqrt(re+1i*im);	
	else
		[rr, ri]  = root_complex(re,im);
		% wave number
		k = rr+1i*ri;
	end
	% cases are separated by
	% az1 = 3/4*pi*h0^2*omega*(g/h0)^(1/2))/(cd*g)
	%     = 3*pi/4*h0^2/cd*k

	% asymptotic linearisation for limit az1 -> 0
	r_low  = 1/2*cd/h0.^2.*8/(3*Pi).*az1;

	% asymptotic linearisation for limit az1 -> infty
	r_high = sqrt( (4/(3*Pi))*cd.*omega.*qt./(g*h0.^3) );
end

