% 2017-04-03 19:11:48.549446883 +0200
%% Gamma parameter for tidal propagation
%% c.f. Cai 2014
function Gamma = Gamma(obj,k,zeta,lambda,mu,phi,mode)
	if (nargin()<7)
		mode = obj.mode;
	end

	psi   = phi/(mu*lambda);
	switch (mode)
	case {'lorentz'}
		% cai table 3
		% eq 19
		L = friction_coefficient_lorentz(phi);
		Gamma = L(2)/2 - k*zeta*L(1)/(3*mu*lambda);
	case {'dronkers'}
		p     = friction_coefficient_dronkers(-phi);

		% NOTE eq 7 in cai 2016 is simpler
		Gamma = 1/pi*(	-   p(1)*(4/3*k/mu*zeta/lambda) ...
				+   p(2)*(1 + 4/3*k*zeta*psi) ...
				- 2*p(3)*phi*(1 + 2/3*k*zeta*(1/psi + psi)) ...
				+   p(4)*phi^2*(3 + 1/psi^2 + 4*k*zeta*(1/psi - psi/3)) ...
			     );
	case {'godin'}
		% coefficients given in the supplement of cai 2014
		G     = friction_coefficient_godin(phi);
		Gamma =   G(1) ...
		        + G(2)*(mu*lambda)^2 ...
			+ k*zeta*(G(3)*mu*lambda + G(4)/(mu*lambda));

		% without river discharge (B2 in cai 2012 sea level)
		%Gamma = 16/(15*pi) + 32/(15*pi)*mu^2*lambda^2;
	case {'savenije'}
		% Note eq 19-20 in cai 2012 phi is not devided by mu*lambda (ps == phi)!!!
		%psi = phi;
		if (psi < 1) % <=> phi < mu*lambda
			% eq 14 in cai 2014
			Gamma = mu*lambda*(1 + 8/3*zeta*psi + psi^2);
		else
			% eq 5 in cai 2014
			Gamma = mu*lambda*(4/3*zeta + 2*psi + 4/3*zeta*psi^2);
		end
	case {'hybrid'}
		GL = obj.Gamma(obj,k,zeta,lambda,mu,phi,'Lorentz');
		GS = obj.Gamma(obj,k,zeta,lambda,mu,phi,'Savenije');
		Gamma = 1/3*GL + 2/3*GS;
	otherwise
		error('unknown friction mode');
	end % switch
end % function Gamma

