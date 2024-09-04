% Sat  1 May 18:37:47 CEST 2021
% TODO implement for higher components
% c.f. lamb-1916, §174
function E = energy(obj,h0,w0,z1,u1)
	g = obj.g;
	omega_1 = obj.omega(1);
	% TODO add 35 kg more for salt water
	rho_w = obj.rho_w;
	% E = w_0 \int_0^T (g h_0^{1/2} |eta|^2 + h_0^{3/2} |u|^2
	%   =  int_0^T a^2 sin^2 dt = 1/2 a^2
	% 1/4 rho w T c (g |z1|^2 + h_0 |u1|^2)
	T = 2*pi./omega_1;
	% celerity
	c0 = sqrt(g*h0);
	% tidal wave length : lambda = c/f = c*T
	lambda = c0*T;
	% E = 0.25*rho_w*T*w0.*sqrt(g*h0).*(g*abs(z1).^2 + h0*abs(u1).^2);
	% E = E_pot + E_kin = 1/2 rho int_L (g eta^2 + h u^2) dx
	% TODO coefficient in lamb is 1/2 not 1/4
	E = 0.25*rho_w*w0.*L.*(g*abs(z1).^2 + h0*abs(u1).^2);
end
