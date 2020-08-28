% 2017-04-04 18:11:50.871153511 +0200
% Karl Kastner, Berlin
%
%% compute friction with the method of Godin
function [uau, obj] = friction_godin(obj,u,Umax)
	if (nargin()<2 || isempty(Umax))
		U   = max(abs(u));
	end

	% chebycheff coeffcients for zero river flow (albeit applied by godin to cases with river flow)
	% c.f. godin 1990, table 1
	% Note: the coefficients do indeed not sum up to 1
	% Note: Godin tries several slightly different sets of coefficients,
	%       of which the Chebysheff set is best
	if (~issym(u))
		% Godins coefficient are just dronkers for zero river flow
		c = obj.friction_coefficient_dronkers(0);
		c = c/pi;
		% c   = [0.3395 0.6791];
	%switch (order)
	%case {3}
	%	c   = [0 0.3395 0 0.6791];
	%case {5}
	%	% 5-term (3 odd)
	%	c   = [0 0.2183  0 1.1641 0 -0.3880];
	%end
	else
		syms A B C;
		c = [A, B, C];
	end

	% chebycheff coeffcients
	% c.f. godin 1990, table 1
	% Note: the coefficients do indeed not sum up to 1
	u_  = u/U;
	uau = 0;
	% even coefficients are zero
	for idx=2:2:length(c)
		%uau = U^2*(c(1)*u_ + c(2)*u_.^3);
		uau = uau + c(idx)*U^2*u_.^(idx-1);
	end
end % River_Tide/friction_godin

