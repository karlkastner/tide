% Thu  7 Jul 16:38:20 CEST 2016
% Karl Kastner, Berlin
%
%% damping modulus of the river tide
%% c.f. Jay and Kukula
%
%% function r = damping_modulus(obj,h0,b,Qr)
%  old: function r = damping_modulus(obj,cD,h0,b,Qr,T)
function [r,obj] = damping_modulus(obj,h0,b,Qr)
	g     = obj.g;
	% Constatn.gravity;
	omega = obj.omega;
	cD    = obj.cD;
	% eq. 7 K&J 2003a and eq. 3 K&J 2003b
	% r = - kappa
	r     = -sqrt(cD.*Qr.*omega / (2*g*h0.^3.*b));
end % River_Tide_JK/damping_modulus

