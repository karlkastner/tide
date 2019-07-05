% Mon  4 Mar 18:47:29 CET 2019
%% discharge amplitude
function [x,qmag,phi,s,c,sl,sr,cl,cr] = discharge_amplitude(obj,id)
		s0 = obj.s(id,:);
		c0 = obj.c(id,:);
		L  = obj.L(id);
		k  = obj.k(id);
		r  = obj.r(id);

		% grid for evalutation
		nx = max(2,L/obj.dx);
		x  = linspace(0,L,nx);

		% sin(o*t - kx) = sin(ot)cos(kx) + cos(ot)*sin(-kx)
		% cos(o*t - kx) = sin(ot)sin(kx) + cos(ot)*cos( kx)

		% sin(ot) of left going wave
		sl = (   s0(1)*exp(-r*x).*cos( k*x) ...
		       + c0(1)*exp(-r*x).*sin(+k*x) );
		% sin(ot) of right going wave
		sr = (   s0(2)*exp(-r*(L-x)).*cos( k*(L-x)) ...
		       + c0(2)*exp(-r*(L-x)).*sin(+k*(L-x)) );
		% sin(ot) of combined wave
		s  = sl+sr;
		
		cl = (   s0(1)*exp(-r*x).*sin(-k*x) ...
		       + c0(1)*exp(-r*x).*cos( k*x) );
		cr = (   s0(2)*exp(-r*(L-x)).*sin(-k*(L-x)) ...
		       + c0(2)*exp(-r*(L-x)).*cos( k*(L-x)) );
		c  = cl+cr;

		% magnitude and phase angle
		qmag = hypot(s,c);
		phi  = atan2(s,c);
end % amplitude

