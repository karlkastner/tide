% Tue  4 Apr 12:12:34 CEST 2017
%% friction coefficient according to Godin
%% these coefficients are identical to Dronker's for U_R = phi = 0
%%
%% function G = friction_coefficient_godin(obj,phi)
function G = friction_coefficient_godin(obj,phi)
	G(1) = 16/(15*pi)*(1 + 2*phi + 7*phi^2)/(1+phi);
	G(2) = 32/(15*pi)*(1/(1+phi));
	G(3) = 128/(15*pi)*phi/(1+phi);
	G(4) = 64/(15*pi)*(1/3*phi + 2/3*phi^2 + phi^3)/(1+phi);
end % River_Tide/friction_coefficient_godin

