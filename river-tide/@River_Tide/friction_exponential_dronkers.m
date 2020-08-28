% Thu 23 Mar 10:51:00 CET 2017
% Karl Kastner, Berlin
%
%% friction coefficicients for the frequency components computed by Dronkers method
%% c.f. Dronker's 1964 eq 8.2 and 8.4
%% Note: Cai dennominates alpha as phi
%%
%% function [c uau uau_ p] = friction_trigonometric_dronkers(u,dp,Umid,Uhr,order,psym)
function [c, uau, uau_, p] = friction_exponential_dronkers(obj,u,dp,Umid,Uhr) %,order,psym)
	if (nargin() < 4 || isempty(Umid))
		[Urange, Umid] = fourier_range(cvec(u),[0;dp]);
		Uhr = 1/2*Urange;
	end
	if (~obj.issym)
		pi_ = pi;
	else
		syms pi_;
	end

	% Note that here midrange and 1/2 range are required,
	% not amplitude and mean
	alpha = Umid/Uhr;
	p     = -obj.friction_coefficient_dronkers(alpha);

	[up] = fourier_power(u,dp);
	
	% for mean, diurnal, semidiurnal, terdiurnal
	p = p/pi_;
	for idx=1:4
		uau(idx).c = 0;
		for jdx=1:obj.opt.friction_order
			uau(idx).c = uau(idx).c + (p(jdx+1)*Uhr.^(2-jdx)*up(idx,1,jdx));
		end
	end % for idx
	% add constant term for mean
	uau(1).c = uau(1).c + (p(1)*Uhr.^2);
	uau_ = [];

	c = [uau(1).c, uau(2).c, uau(3).c, uau(4).c].';
end % River_Tide/friction_exponential_dronkers

