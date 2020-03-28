% ode-coefficients of the SWE jacobian transformed to the wave equation
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
end

