% 2020-03-15 18:31:25.716017562 +0800
% Karl Kastner, Berlin
%
%% Jacobian matrix indices of the Shallow-Water Equation
%% d(A,Q)/dt + J(A,Q) d(A,Q)/dx = forcing terms
%%
%%  [       0,      1][dA/dx] 
%%  [-Q^2/A^2, 2 Q/A ][dQ/dx]
%
function f = odefun_swe_jacobian(obj,f,k,h0,w0)
	g = obj.g;
	omega = obj.omega;

	% ode = - 1/(g h w) d^2/dt^2 Q + 1/w d^2/dx^2 Q = 0

	% Q''
	f(:,1) = f(:,1) + 1./w0;

	% Q' 
	% 0

	% Q
	f(:,3) = f(:,3) - (1i*k*omega).^2./(g.*h0.*w0);

	% constant
	% 0
end % odefun_swe_jacobian

