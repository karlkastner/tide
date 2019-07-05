% Wed  6 Jul 13:39:00 CEST 2016
% Karl Kastner
%% tidal range
%%
%% based on Horrevoets/Savenije, 2004
%%
%% H0   : tidal range at river mouth
%% h0   : initial water depth
%% v    : velocity scale
%% b    : convergence length
%% sine : phase lag
%% K    : Mannings coefficient
%% Q_r   : river discharge

function [x eta eta_] = savenije_tidal_range1(X,T,eta0,h0_fun,b0,L_b,K,U_t,sine,I,Q_r)
	y0     = 1;
	H0     = 2*eta0;
	[x yH] = ode23s(@(x,y) ydot(x,y,T,H0,h0_fun,b0,L_b,K,U_t,sine,I,Q_r),X,[y0 H0]);
	eta    = eta0*yH(:,1);
	eta_   = 0.5*yH(:,2);
end

function [yHdot] = ydot(x,yH,T,H0,h0_fun,b0,L_b,K,U_t,sine,I,Q_r)
	H = yH(2);
	y = yH(1);

	% acceleration by gravity
	g  = 9.81;
	
	% depth
	h0 = h0_fun(x);
	%H  = y*H0;

	% cross section area
	% note: horrevoets denotes it with A, but it is A0	
	A  = h0.*b0.*exp(-x/L_b);
	b = L_b;

	% cross sectional average velocity
	U_r = Q_r/A;

	a = L_b;
	delta = 0;
	eps = savenije_phase_lag(T,h0,U_t,U_r,a,delta);
	sine = sin(eps);

	% friction factor, eq 14
	f  = 1./(1 - ((1.33*H)/(2*h0)).^2);
	fp = g/(K^2*h0^(1/3))*f;


	% wave celerity
	% Typo in savenije, c_2 should be c^2
	c0 = sqrt(g*h0);

	phi = U_r/U_t;

	% correction factor of wave celertity, eq. 15
	% note, this is theta in Savenije 2012
	gamma = 1 - (sqrt(1+H/(2*h0))-1)*phi/sine;

	% note: condition differs to Savenije
	if (U_t*sin(eps) > U_r)
		% equation 12
		Hdot = 1./(g/(2*c0*U_t*sine) - phi/(sine*2*h0) + gamma/H) ...
			.* (  gamma/b ...
                              - fp*U_t*sine/(h0*c0)*(1 + 4/3*H/h0*phi/sine ...
                                                         + phi^2/sine^2) ...
			- I/h0);
	else
		% equation 13
		Hdot = 1./(g/(2*c0*U_t*sine) - phi/(sine*2*h0) + gamma/H) ...
			.* ( gamma/b ...
                             - fp*U_t*sine/(h0*c0)*(2/3*H/h0 + 2*phi/sine ...
						+ 2/3*H/h0*phi^2/sine^2 ) ...
			   - I/h0);
	end

	% Tidal Froude number
	a  = U_t*sine/c0*2*h0/H0;
	at = a./(1 - c0.*U_r./(g*h0));
	
	% eq. 18
	beta1 = 1./( gamma./L_b ...
		  - fp*U_t*sine/(h0*c0) ...
		  .*(1 + 4/3*y.*H0./h0.*phi./sine ...
		       + phi^2/sine.^2));
	% eq 20a, note Savenije has I/h instead of I/h0
	ydot = 1./(1./at + gamma./y).*(1./beta1 - I./h0);

	yHdot = [ydot; Hdot];
end % function ydot

