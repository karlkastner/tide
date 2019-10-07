% Tue  4 Apr 15:55:45 CEST 2017
%% friction computed by the method of Lorent'z
%% expressed as coefficients of the frequency components (trigonometric form)
%% 
% c.f. dronkers 4.23
function [c, uau] = friction_trigonometric_lorentz(obj,u,dp)
	% phase
	phi = [0, 0, dp];

	u0 = u(1);

	% note that dronkers defines U1 = 1/2*u1
	u1 = u(2)*cos(phi(2));

	phi   = u0/u1;
	if (abs(phi)<1)
		% river to tidal ratio
		gamma = acos(phi);
		k00 = (2+cos(2*gamma))*(2-4*gamma/pi)+6/pi*sin(2*gamma);
		k10 = 6/pi*sin(gamma) + 2/(3*pi)*sin(3*gamma) + (4-8/pi*gamma)*cos(gamma);

		% constant term of the friction
		% 4.29 in dronkers
		C0  = k00*1/4*abs(u1)^2;

		% leading oscillating term (diurnal or semidiurnal, depends on user definition of u)
		C1  = k10*1/4*abs(u1);
	else
		% note: the computation with k is not save, because it leads to divisiob by zero if 0==u1
		% 4.34 in dronkers
		C0 = -(-u0^2 - 1/2*abs(u1)^2);
		C1 = u0;
	end

	% constant term
	uau(1).a   = C0;

	% diurnal term
	uau(2).c = 2*C1*u1;
	uau(2).a = real(uau(2).c);
	uau(2).b = imag(uau(2).c);
		
	uau(3).c = NaN;
	uau(3).a = NaN;
	uau(3).b = NaN;

	uau(4).c = NaN;
	uau(4).a = NaN;
	uau(4).b = NaN;
	
	c=[uau(1).a uau(2).a uau(2).b uau(3).a uau(3).b uau(4).a uau(4).b];

	% semidiurnal term, complicated, see dronkers
end % River_Tide/friction_trigonometric_lorentz

