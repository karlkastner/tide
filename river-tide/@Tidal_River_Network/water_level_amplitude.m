% Fri 22 Feb 14:22:25 CET 2019
%% predict the surface elevation amplitude
function [x,zmag,phi,s,c,sl,sr,cl,cr] = water_level_amplitude(obj,id)
		s0 = obj.z1s(id,:);
		c0 = obj.z1c(id,:);
		L  = obj.L(id);
		k  = obj.k(id);
		r  = obj.r(id);
		w  = obj.width(id);
		o  = obj.omega(1);

		% grid for evalutation
		nx = max(2,L/obj.dx);
		x  = linspace(0,L,nx);

		% sin(o*t - kx) = sin(ot)cos(kx) + cos(ot)*sin(-kx)
		% cos(o*t - kx) = sin(ot)sin(kx) + cos(ot)*cos( kx)

		% sin(ot) of left going wave
		% z = -1/iow dQ/dx
		sl = -1/(1i*o*w)*exp(-r*x).*(  s0(1)*(-r*cos(k*x) - k*sin(k*x)) ...
			                    + c0(1)*(-r*sin(k*x) + k*cos(k*x)) );


		% sin(ot) of right going wave
		sr = -1/(1i*o*w)*exp(-r*(L-x)).*(  s0(2)*(r*cos( k*(L-x)) + k*sin(k*(L-x))) ...
			  			+ c0(2)*(r*sin(+k*(L-x)) - k*cos(k*(L-x))));

		% sin(ot) of combined wave
		s  = sl+sr;
		
		cl = -1/(1i*o*w)*exp(-r*x).*(   s0(1)*(-r*sin(-k*x) -k*cos(k*x)) ...
		                  + c0(1)*(-r*cos( k*x) -k*sin(k*x)));

		cr = -1/(1i*o*w)*exp(-r*(L-x)).*(   s0(2)*(+r*sin(-k*(L-x)) + k*cos(k*(L-x))) ...
		                      + c0(2)*(+r*cos( k*(L-x)) + k*sin(k*(L-x))));

		% cos(ot) of combined wave
		c  = cl+cr;

		% magnitude and phase angle
		zmag = hypot(s,c);
		phi  = atan2(s,c);
end % water_level_amplitude

