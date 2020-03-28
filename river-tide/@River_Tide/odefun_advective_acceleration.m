% Sun  8 Oct 13:08:39 CEST 2017
% advective acceleration, ignoring the 1/h nonlinearity
%
% aa = d/dx (Q^2/A)
%
function f = odefun_advective_acceleration(obj,f,k,Q,QQ,h0,dh0_dx,w0,dw0_dx,D1_dx)
	g = obj.g;
	omega  = obj.omega;
	flag   = obj.flag;

	s      = -1./(g.*h0.*w0);
	dA0_dx = h0.*dw0_dx + w0.*dh0_dx;
	A0     = h0.*w0;

	% self-interaction, first-derivative, (d/dx(Q0*Qk) = (dQ0/dx*Qk + Q0*dQk/dx)
	f(:,2) = f(:,2) + s.*(2./A0.*Q(:,1));
	
	% self-interaction, zero-derivative
	f(:,3) = f(:,2) + s.*(-2.*Q(:,1)./(A0.*A0).*dA0_dx + 2./A0.*(D1_dx.*Q(:,1)));

	% non-linear interaction, forcing
	f(:,4) = f(:,4) + s.*(-QQ(:,k+1)./(A0.*A0).*dA0_dx + 1./A0.*(D1_dx.*QQ(:,k+1)));
end

