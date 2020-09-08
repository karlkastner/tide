% 2017-04-04 18:11:50.871153511 +0200
% Karl Kastner, Berlin
%
% friction coefficient according to lorentz
% only applicable to the case of no river flow
function [uau, obj] = friction_lorentz(obj,u,Umax)
	if (nargin()<2 || isempty(Umax))
		U   = max(abs(u));
	end

	if (~issym(u))
		% Godins coefficient are just dronkers for zero river flow
		c = obj.friction_coefficient_lorentz();
		c = c/pi;
	else
		syms A B C;
		c = [A, B, C];
	end

	u_  = u./U;
	uau = 0;
	% even coefficients are zero
	for idx=2:2:length(c)
		%uau = U^2*(c(1)*u_ + c(2)*u_.^3);
		uau = uau + c(idx)*U^2*u_.^(idx-1);
	end
end % River_Tide/friction_lorentz

