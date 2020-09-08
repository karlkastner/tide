% Wed  5 Apr 10:20:14 CEST 2017
% Karl Kastner, Berlin
%
%% friction computed by the method of Godin
%% expressed as coefficients of the frequency components (trigonometric form)
%%
%% Chebycheff coeffcients for zero river flow
%% (albeit applied by Godin to cases with river flow)
%% c.f. godin 1990, table 1, column Ch
%% Note: the coefficients do indeed not (exactly) sum up to 1
%% Note: Godin tries several slightly different sets of coefficients,
%%       of which the Chebysheff set is best
%%
%   function [c, uau] = friction_trigonometric_godin(obj,u,dp,Umax)
function [c, uau] = friction_trigonometric_godin(obj,u,dp,Umax)
	if (nargin() < 3 || isempty(Umax))
		[Urange, Umid] = fourier_range(cvec(u),[0;dp]);
		Umax = 1/2*Urange + abs(Umid);
	end
	
	pi_ = obj.pi;
	if (1) %~issym(u))
		% Godins coefficient are just dronkers for zero river flow
		%c3   = [0, 0.3395,  0, 0.6791, 0,  0     ];
		%c5   = [0, 0.2183,  0, 1.1641, 0, -0.3880];
		p = obj.friction_coefficient_dronkers(0);
	else
		syms A B C D E
		p = [A B C D E];
	end

	% was devided by Umax and below in loop power was not adjusted
	[up] = fourier_power(u,dp);

	% for mean, diurnal, semidiurnal, terdiurnal
	p = p/pi_;
	for idx=1:4
		uau(idx).a = 0;
		uau(idx).b = 0;
		for jdx=1:obj.opt.friction_order
			% TODO, power was just 2, not 2-jdx
			uau(idx).a = uau(idx).a + (p(jdx+1)*Umax^(2-jdx)*up(idx,1,jdx));
			uau(idx).b = uau(idx).b + (p(jdx+1)*Umax^(2-jdx)*up(idx,2,jdx));
		end % for jdx
	end % for idx
	% add constant term for mean
	uau(1).a = uau(1).a + (p(1)*Umax.^2);

	c = [uau(1).a, uau(2).a, uau(2).b, uau(3).a, uau(3).b, uau(4).a, uau(4).b].';
end % River_Tide/friction_trigonometric_godin

