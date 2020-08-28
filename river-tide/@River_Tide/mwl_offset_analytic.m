% Wed 18 Apr 15:26:01 CEST 2018
% analytic solution:
%	    prasmatic channel with uniform flow (tidal average)
%           z0' << h0
% 	    Q0 >> Qt (z1 << h0)
function z0t = mwl_offset_analytic(obj,x,z10,h0,w0,Cd,Q0)
	g  = Constant.gravity;
%	x  = obj.x(cdx);
%	w0 = obj.width(cdx,x(1));
%	h0 = obj.h0(cdx,x(1));
%	Q0 = obj.Q(0,cdx,x(1));
%	Q0 = obj.Q0_;
	omega = obj.omega;
%	Cd = obj.cd(cdx,x(1));
	r1  = sqrt(omega*Cd*Q0/(g*w0*h0.^3));
	k   = (1-1i)*r1;
	%k = r1*sqrt(2);
	Q10 = abs(1i*omega*w0/k)*z10;
	% r0  = 1/2*(3*Q0^2*Cd)./(g*h0^4*w0^2); <- from script
	% thesis seems factor 3 too small
	r0 = 1/2*Cd*Q0^2/(g*w0^2*h0^4); % <- from thesis
	% correction:
	%r0 = 3*r0;
	%z0t =     Q10^2*cd*h0/(2*(3*cd*Q0^2 - 2*g*r1*h0^4*w0^2)) ...
	% from thesis
	z0t =     (Q10/Q0)^2*h0/6*r0/(r0-r1) ...
		  .*(exp(-2*r1*x) - exp(-2*r0*x));
end

