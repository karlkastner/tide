% Thu 10 Oct 11:36:29 PST 2019
function [rmse,res] = check_momentum(obj)
	% TODO use flags for 1/h and advective acceleration
	x   = obj.x;
	z1  = obj.z_(:,2);
	Q0  = obj.Q_(:,1);
	Q1  = obj.Q_(:,2);
	w   = obj.width(x);
	h0  = obj.h0(x);
	cD  = obj.cd(x);
	dz1_dx = derivative1(x,z1);
	aa  = 0;
	aQQ = 0;
	n   = 100;
	t   = (0:n-1)/n*T;
	% time series of discharge
	Q_   = Q0*ones(1,n) + Q1*exp(1i*omega*t);
	aQQ  = abs(Q_).*Q_;
	% omega_1 frequency component
	aQQ1 = aQQ*exp(1i*omega*t)';
	res  = 1i*omega*Q1 ...  % acceleration in time
		+ aa ...        % advective acceleration
		+ g*dz1_dx ...  % pressure gradient
		+ cD*aQQ/(w^2.*h0.^3); % friction
	% source terms for width and bed level?
end % check_momentum

