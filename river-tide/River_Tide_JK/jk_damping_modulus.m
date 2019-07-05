% Thu  7 Jul 16:38:20 CEST 2016
% Karl Kastner, Berlin
%
%% damping modulus of the river tide
%% c.f. Jay and Kukula
%
function r = jk_damping_modulus(cD,h0,b,Qr,T)
	g     = Constatn.gravity;
	omega = 2*pi./T;
	% eq. 7 K&J 2003a and eq. 3 K&J 2003b
	% r = - kappa
	r     = -sqrt(cD.*Qr.*omega / (2*g*h0.^3.*b));
end

