% Tue  4 Apr 11:18:46 CEST 2017
%% Loren't friction coefficient
%% c.f. cai
%% c.f. dronkers 
function L = friction_coefficient_lorentz(obj,phi)
	if (phi < 1)
		% alpha is gamma in dronkers
		alpha = acos(-phi);
		% eq 17 in cai, 4. 30 in dronkers
		L(1) = (2 + cos(2*alpha))*(2 - 4*alpha/pi) + 6/pi*sin(2*alpha);
		% eq 18 in cai, 4.32 in dronkers
		L(2) = 6/pi*sin(alpha) + 2/(3*pi)*sin(3*alpha) + (4-8/pi*alpha)*cos(alpha);
	else
		% eq 20
		L(1) = -2 - 4*phi^2;
		L(2) =      4*phi;
	end % else of if
end % fcl

