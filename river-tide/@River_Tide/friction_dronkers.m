% Thu 23 Mar 09:51:57 CET 2017
% Karl Kastner, Berlin
%
%% friction determined by Dronker's method
%%
%% input :
%%         u    : velocity time series
%%	   Umid : arithmetic mean of mininmum and maximum velocity
%%                (not the mean of the velocity, usually non-zero even without river flow)
%%         Uhr  : half-range of the velocity, less than the sum of
%%                the frequency amplitudes, except at perigean spring tides
%%
%% function [uau_sum uau p] = friction_dronkers(u,Umid,Uhr,order)
function [uau_sum, uau, p, obj] = friction_dronkers(obj,u,Umid,Uhr)
	if (nargin() < 2||isempty(Uhr))
		Umin = min(u);
		Umax = max(u);
		Uhr  = 1/2*(Umax - Umin);
		Umid = 1/2*(Umax + Umin);
	end
%	if (nargin()<4 || isempty(order))
%		order = 3;
%	end

	% note: in dronkers, U and not Q is used
	% Ut (Qa) : one half of tidal velocity range
	% Ur (Qb) : river velocity

	% dronkers 8.2 and 8.4 (dronkers uses alpha : phi)

	if (issym(u)) %obj.issymbolic) %nargin()>4 && psym)
		syms p0 p1 p2 p3
		p = [p0, p1, p2, p3];
		syms pi_
	else
		alpha = Umid/Uhr;
		p     = -obj.friction_coefficient_dronkers(alpha);
		pi_ = pi;
	end

	order = obj.opt.friction_order;
	for idx=1:order+1
		% uau = 1/pi (p0 U^2 + p1 u U + p2 u^2/U + p3 u^3/U^2)
		% uau = 1/pi (p0 U + p1 u + p2 u^2 + p3 u^3/U)
		uau(:,idx) = 1/pi_*p(idx)*Uhr.^(3-idx)*u.^(idx-1);
	end
	uau_sum = sum(uau,2);
end % River_Tide/friction_dronkers

