% Fri 27 Jan 12:38:16 CET 2017
%% energy transport of a tidal wave
function [E, Et] = energy_transport_1d(t,h,Q,u,zs)
	% TODO no magic numbers
	rho = 1000;

	dt = t(2)-t(1);
	n  = round(25/24*1/dt);

	% daily average depth
	% should hbar here better computed with respect to mean sea level?
	hbar = meanfilt1(h',n)';
	Qbar = meanfilt1(Q',n)';
	ubar = meanfilt1(u',n)';

	% surface level variation
	ht = h - hbar;
	ut = u - ubar;
	Qt = Q - Qbar;

	% total energy transport
%	E = rho*Q.*(1/2*u.^2 + ht);
	E = rho*Q.*(1/2*u.^2 + zs);

	% tidal transport, without river flow
	Et = rho*Qt.*(1/2*ut.^2 + ht);
end % energy transport

