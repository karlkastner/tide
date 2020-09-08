% Tue  4 Apr 11:18:46 CEST 2017
% Karl Kastner, Berlin
%
%% coefficients of the Fourier expansion of the signed square of the |Q|Q
%% of the friction term
%%
%% Lorent'z used this first for the case of no river flow
%%
%% identical to Dronker's coefficient for zero river flow
%% and a single frequency component
%% c.f. Cai
%% c.f. Dronkers (gamma = alpha) 
%%
%% note difference in coefficients due to different definitions:
%% definition used here:
%%	Q = Q0 + 1/2*(sum_k Q_k e(k iwt) + conj(Q_k e(k iwt)))
%% but Dronkers defines
%%	Q = Q + sum_k Q_k e(k iwt)
%%
%% function L = friction_coefficient_lorentz(obj,phi)
function L = friction_coefficient_lorentz(obj,phi)
	pi_ = obj.pi;
	if (phi < 1)
		% Drokers defines alpha in eq 4.22 as the moments of slack water: 
		% alpha = acos(Q0/(2*abs(Q1))

		% alpha is gamma in Dronkers
		alpha = acos(-phi);

		% eq 17 in cai, p 274, 4.30 in Dronkers 1964
		L(1) = (   (2 + cos(2*alpha)).*(2 - 4*alpha/pi_) ...
                         + 6/pi_*sin(2*alpha) ...
		       );

		% eq 18 in cai, 4.32 in dronkers
		L(2) = (    6/pi_*sin(alpha) ...
			+ 2/(3*pi_)*sin(3*alpha) ...
			+ (4-8/pi_*alpha)*cos(alpha) ...
		       );
	else
		% eq 20, 4.35 in Drokers, 
		% note that dronkers has factor 2, and second sign is inverted
		L(1) = -2 - 4*phi^2;
		L(2) =      4*phi;
	end % else of if
end % River_Tide/friction_coefficient_Lorentz

