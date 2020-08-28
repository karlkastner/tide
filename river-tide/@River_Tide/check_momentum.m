% Thu 10 Oct 11:36:29 PST 2019
function [rmse,res] = check_momentum(obj,obj,x,w,h0,Q0,zt,Qt,omega)	
	% TODO use flags for 1/h and advective acceleration
	aa  = 0;
	aQQ = 0;
	n   = 100;
	t   = (0:n-1)/n*T;
	% time series of discharge
	Q_   = Q0*ones(1,n) + Qt*exp(1i*omega*t);
	aQQ  = abs(Q_).*Q_;
	dz1_dx = derivative1(x,zt);
	% omega_1 frequency component
	aQQt = aQQ*exp(1i*omega*t)';
	res  = 1i*omega*Qt ...  % acceleration in time
		+ aa ...        % advective acceleration
		+ g*dz1_dx ...  % pressure gradient
		+ cD*aQQ/(w^2.*h0.^3); % friction
	% source terms for width and bed level?
end % check_momentum

