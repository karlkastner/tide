% 2020-03-15 18:31:25.716017562 +0800
% Karl Kastner, Berlin
%
%% Jacobian matrix of the Shallow-Water Equation
%% d(A,Q)/dt + J(A,Q) d(A,Q)/dx = forcing terms
%%
%% dA/dt  + [       0,      1][dA/dx] = [f_c]	(c)
%% dQ/dt    [-Q^2/A^2, 2 Q/A ][dQ/dx]   [f_m]	(m)
%%
%% dm/dt 
%% d^2Q/dt^2  - d/dt Q^2/A^2 - 2 d/dt (Q/A) dQ/dx = d/dt f_m
%% d^2Q/dt^2  - d/dt Q^2/A^2 - 2 (1/A dQ/dt - Q/A^2 dA/dt) dQ/dx = d/dt f_m
%% d^2Q/dt^2  - d/dt Q^2/A^2 - 2 (1/A dQ/dt + Q/A^2 dQ/dx) dQ/dx = d/dt f_m
%%
%% ode = - 1/(g h w) d^2/dt^2 Q + 1/w d^2/dx^2 Q = 0
%%     = -1 (iko)^2/(g h w) Q + 1/w Q''
%%     =    (k^2 o^2)/(g h w) Q + 1/w Q''
function f = odefun_swe_jacobian(obj,f,k,h0,w0)
	g     = obj.g;
	omega = obj.omega;


	% Q''
	f(:,1) = f(:,1) + 1./w0;

	% Q' 
	% f(2) = 0

	% Q
	f(:,3) = f(:,3) + (k*k*omega*omega)./(g.*h0.*w0);
	%f(:,3) = f(:,3) - (1i*k*omega).^2./(g.*h0.*w0);

	% constant
	% f(4) = 0
end % odefun_swe_jacobian

