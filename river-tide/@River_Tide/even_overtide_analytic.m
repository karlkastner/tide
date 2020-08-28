% Wed 25 Apr 15:09:11 CEST 2018
function z2 = even_overtide_analytic(obj,x,z10,h0,w0,Q0,Cd,omega)
	g     = Constant.gravity;	
%	x  = obj.x(cdx);
%	cD = obj.cd(cdx,x(1));
%	w0 = obj.width(cdx,x(1));
%	h0 = obj.h0(cdx,x(1));
%	Q0 = obj.Q0(cdx,x(1));
%	omega = obj.omega;

	r1 = sqrt(omega*Cd*Q0/(g*w0*h0^3));
	k1 = (1-1i)*r1;
	k2 = sqrt(2)*k1;
	r2 = sqrt(2)*r1;
	z2 = omega*w0/(8*r1*Q0)*z10^2*(exp(-2i*k1*x) - exp(-1i*k2*x));
end

