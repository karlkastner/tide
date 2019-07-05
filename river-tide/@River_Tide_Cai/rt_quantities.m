% 2017-04-03 19:16:46.345006529 +020
%% determine the quantities that determine the tidal propagation
%% c.f. Cai
%% 
%% Note: this computes 4 unknowns following Cai, however,
%% lambda, mu and epsilon can be substituted
%% making it an equation in one unknown (delta) only
function [delta, lambda, mu, e] = rt_quantities(obj,zeta,gamma,chi,rs,phi)

	% k=1 : consider depth variation in friction term
	k = 1;

	% start value of iteration
	% attenuation
	% note : if start value > -1, then the iteration does not converge to the correct value for large phi!
	delta  = -1;

	% iterate until convergence
	% TODO use pade
	iter = 0;
	while (1)
		iter = iter+1;
	
		
		% (inverse) celerity scale
		% tab. 1, c = c0/lambda;
		% cai 2016, eq 15
		lambda = sqrt(1-delta*(gamma-delta));

		% phase lag
		% cai 2016: eq. 16
		e      = atan(lambda/(gamma-delta))

		% velocity number
		% tab 1: mu = rs*zeta*c0;
		% cai 2016, eq. 14
		mu     = sin(e)/lambda;

		if (1)
		% correction factor of the wave celerity
		% eq. 13 in cai 2014
		% tab 1 in cai 2016
		% note: cai 2012 20-21 misses both 1/(mu*lambda) in theta and beta terms
		theta  = 1 - (sqrt(1+zeta)-1)*phi/(mu*lambda);

		% correction factor of the tidal froude number
		% tab 1 in cai 2016
		% eq 12 in cai 2014
		beta   = theta - rs*zeta*phi/(mu*lambda);
		else
			theta  = 1 - (sqrt(1+zeta)-1)*phi;
			beta   = theta - rs*zeta*phi;
		end
		%theta=1;
		%beta=1;

		% friction
		Gamma = obj.Gamma(k,zeta,lambda,mu,phi);

		% Damping number
		% eq. 11 in cai 2014
		% eq. 6 in Cai 2016
		delta_new = mu^2/(1+mu^2*beta)*(gamma*theta - chi*mu*lambda*Gamma);

		% over relaxation to make iteration converge
		delta_old = delta;
		delta = obj.relax*delta_new + (1-obj.relax)*delta;

		delta_new-delta_old
		if (abs(delta_new-delta_old) < obj.abstol ...
			|| abs(delta_new-delta_old) < obj.reltol*abs(delta_old))
			break;
		end % break
		if (~isfinite(delta) || ~isreal(delta))
			warning('non-finite or imaginary value');
			return;
		end
		if (iter > obj.maxiter)
			error('maximum number of iterations reached');
		end

	end %while
end % rt_quantities

