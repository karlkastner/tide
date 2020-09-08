% Thu 23 Mar 10:51:00 CET 2017
% Karl Kastner, Berlin
%
%% friction computed by the method of Dronkers
%% expressed as coefficients for the frequency components
%% c.f. dronkers 1964 eq 8.2 and 8.4
%% Note: Cai dennominates alpha as phi
% function [c, uau, uau_, p] = friction_trigonometric_dronkers(u,dp,Umid,Uhr,order,psym)
function [c, uau, uau_, p] = friction_trigonometric_dronkers(obj,u,dp,Umid,Uhr)
	if (nargin() < 4 || isempty(Umid))
		[Urange, Umid] = fourier_range(cvec(u),[0;dp]);
		Uhr = 1/2*Urange;
	end
	pi_ = obj.pi;

	% Note that here midrange and half-range are required,
	% not amplitude and mean
	alpha = Umid/Uhr;
	p     = -obj.friction_coefficient_dronkers(alpha);

	[up] = fourier_power(u,dp);
	
	% for mean, diurnal, semidiurnal, terdiurnal
	p = p/pi_;
	for idx=1:4
		uau(idx).a = 0;
		uau(idx).b = 0;
		for jdx=1:obj.opt.friction_order
			uau(idx).a = uau(idx).a + (p(jdx+1)*Uhr.^(2-jdx)*up(idx,1,jdx));
			uau(idx).b = uau(idx).b + (p(jdx+1)*Uhr.^(2-jdx)*up(idx,2,jdx));
		end % for jdx
	end % for idx
	% add constant term for mean
	uau(1).a = uau(1).a + (p(1)*Uhr.^2);

	c = [uau(1).a, uau(2).a, uau(2).b, uau(3).a, uau(3).b, uau(4).a, uau(4).b].';

	uau_ = [];
end % friction_trigonometric_dronkers

