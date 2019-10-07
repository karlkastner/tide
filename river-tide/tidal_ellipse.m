% Tue 17 Apr 11:30:33 CEST 2018
%% tidal ellipse, numerical ode solution
function [xt]   = tidal_ellipse(x,u0,u1,omega1,x0)
	% normalize angle
	u1_ = interp1(x,u1,x0);
	u1  = u1/exp(1i*(angle(u1_)));

	Ti      = 2*pi/omega1;
	[T, xt] = ode23s(@ufun,[0,Ti],x0);

function u_ = ufun(t_,x_)
	u0_ = interp1(x,u0,x_);
	u1_ = interp1(x,u1,x_);
	u_  = u0_ + real(u1_*exp(1i*omega1*t_));
end

end

