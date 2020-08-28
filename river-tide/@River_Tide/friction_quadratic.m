% Thu 23 Mar 09:51:57 CET 2017
% Karl Kastner, Berlin
%
%% friction determined by Dronker's method
% function [uau_sum uau p] = friction_dronkers(u,Umid,Uhr,order)
function [uau_sum, uau, p, obj] = friction_quadratic(obj,u)
	%order = 2;
	%for idx=1:order+1
		% uau = 1/pi (p0 U^2 + p1 u U + p2 u^2/U + p3 u^3/U^2)
		% uau = 1/pi (p0 U + p1 u + p2 u^2 + p3 u^3/U)
	uau = u.*u;
	%end
	uau_sum = sum(uau,2);
end % friction_dronkers

