% Wed  5 Apr 10:20:14 CEST 2017
% Karl Kastner, Berlin
%
%% friction computed by the method of godin
%% expressed in fourier series coefficients
function [c, uau] = friction_trigonometric_godin(obj,u,dp,Umax)
	if (nargin() < 3 || isempty(Umax))
		[Urange, Umid] = fourier_range(cvec(u),[0;dp]);
		Umax = 1/2*Urange + abs(Umid);
	end
	% chebycheff coeffcients for zero river flow (albeit applied by godin to cases with river flow)
	% c.f. godin 1990, table 1
	% Note: the coefficients do indeed not sum up to 1
	% Note: Godin tries several slightly different sets of coefficients,
	%       of which the Chebysheff set is best
	
	if (obj.issym)
		syms pi_
		syms phi
	else
		phi = 0;
		pi_ = pi;
	end
	if (1) %~issym(u))
		% Godins coefficient are just dronkers for zero river flow
		%c3   = [0, 0.3395,  0, 0.6791, 0,  0     ];
		%c5   = [0, 0.2183,  0, 1.1641, 0, -0.3880];
		phi_ = 0;
		c = obj.friction_coefficient_dronkers(phi_);
		%c = subs(c,phi,0);
	else
		syms A B C D E
		c = [A B C D E];
	end

	[up] = fourier_power(u/Umax,dp);

	c = c/pi_;
	for idx=1:4
		uau(idx).a = 0;
		uau(idx).b = 0;
		for jdx=1:obj.opt.friction_order
			uau(idx).a = uau(idx).a + Umax^2*(c(jdx+1)*up(idx,1,jdx));
			uau(idx).b = uau(idx).b + Umax^2*(c(jdx+1)*up(idx,2,jdx));
		end
	end
	% add constant term for mean
	uau(1).a = uau(1).a + (c(1)*Umax.^2);

	c = [uau(1).a uau(2).a uau(2).b uau(3).a uau(3).b uau(4).a uau(4).b].';
end % friction_trigonometric_godin
