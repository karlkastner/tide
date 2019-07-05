% Tue 28 Nov 14:01:41 CET 2017
%
%% reflection coefficient for gradual varying cross section geometry
%% without damping
%
function [rtotal, r] = reflection_coefficient_gradual(x,h,w,omega)
	g      = 9.81;
	% admittance
	y      = w*sqrt(h);
	c      = sqrt(g*h);
	k      = omega./c;
	% dc_dx = 
	dy_dx  = 1/2*cdiff(y)./cdiff(x);
	% reflection coefficient
	r      = 1./y.*dy_dx.*exp(2i*k.*x);
	dx     = x(2)-x(1);
	rtotal = sum(r*dx);
end
	


