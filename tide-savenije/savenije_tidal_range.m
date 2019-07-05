% Wed 20 Jul 10:36:09 CEST 2016
% Karl Kastner
%% tidal range
%%
%% based on Savenije 2012
%%
%% x    : distance to river mouth
%% eta  : range
%% eta0 : range at river mouth
%% hbar : mean water depth
%% phi  : velocity ratio u_tide/u_river 
%%	 note: this varies in strongly convergent estuaries
%% K    : mannings coefficient
%% I    : residual surface slope I
function [X, eta] = savenije_tidal_range(Xrange,eta0,T,hbar,U_t,U_r,K,I,a,b,delta)
	[X, eta] = ode23s(@(x,eta) deta(x,eta,T,hbar,U_t,U_r,K,I,a,b,delta),Xrange,eta0);
end

function deta = deta(x,eta,T,hbar,U_t,U_r,K,I,a,b,delta)
	% acceleration by gravity
	g = 9.81;

	% eps  : phase lag between hw and hws
	% TODO, sign wrong?
	eps  = savenije_phase_lag(T,hbar,U_t,U_r,a,delta);
	sine = sin(eps);

	phi = U_r/U_t;
	
	% classical wave celerity
	c0 = sqrt(g*hbar);
	
	% Savenije leaves the exact specification of c unclear,
	% Horrevoets does not apply a correction to c, so here it is neither done
	c = c0;

	% tidal froude number (3.23)
	%alpha = mu^2;
	%alpha = sin(eps)^2/(kappa*lambda^2)
	alpha = sin(eps)^2;

	% amplitude to depth ratio
	chi = eta/hbar;

	% adjusted friction factor (3.16)
	f  = 1/(1-(1.33*eta/hbar)^2);
	fp = g/(K^2*hbar^(1/3))*f;

	% correction factor for wave celerity (3.135)
	% note, this is gamma in Horrevoets
	theta = 1 - (sqrt(1+chi)-1)*phi/sine;

	% storage width ratio (assumed unity here, no tidal floodplains) (eq. 2.9)
	r_s = 1;

	% note: condition differs to Horrevoets
	if (eps > phi)
	deta = eta*1./(1/alpha - r_s*phi*chi/sine + theta) ...
		* (   theta/b ...
		    - fp*U_t*sine/(hbar*c)*(1 + 8/3*chi*phi/sine + phi^2/sine^2 ) ...
                    - r_s*I/hbar ...
                  );
	else
	deta = eta*1./(1/alpha - r_s*phi*chi/sine + theta) ...
		* (theta/b ...
		   - fp*U_t*sine/(hbar*c)*(4/3*chi + 2*phi/sine + 4/3*chi*phi^2/sine^2 ) ...
		   - r_s*I/hbar ...
		  );
	end
end

