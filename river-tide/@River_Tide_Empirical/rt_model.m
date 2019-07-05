% Fri 29 Jul 11:58:00 CEST 2016
% Karl Kastner, Berlin
%% select the model for fitting
function obj = rt_model(obj)
	% select function of damping modulus
	switch (obj.model)
	case {1}
		% initial values for nl-regression
		c0 = [0 1e-2];
		fD = @fD1;
	case {2} % excluding self damping
		c0    = [0 1e-2 2];
		fD = @fD2;
	case {3} % including self damping
		c0 = [0 1e-2 2 0.1];
		fD = @fD3;
%	case {4} % including simplified neap-spring-variation
%		c0 = [0 1e-2 2 1e-3];
%		fr = @(c,k,x,x0,Ur,r0,D1) -(   c(1)./(x-x0) ....
%					  + c(2)*Ur.^c(3) ...
%					  + c(4)*r0.^2./sqrt(c(2)*Ur.^c(3)) );
%	case {5} % including separated neap-spring-variation
%		c0 = [0 1e-2 2 1e-3 1];
%		fr = @(c,k,x,x0,Ur,r0,D1)  -(  c(1)./(x-x0) ....
%					     + c(2)*Ur.^c(3) ...
%					     + c(4)*r0.^2.*Ur.^c(5) );
%	case {6} % instationary model w/o neap-spring variation
%		coeff = [0 1e-2 2 0];
%		fr = @(coeff,Ur,r0,dUr_dt) -(coeff(1) + coeff(2)*Ur.^coeff(3) + coeff(4)*dUr_dt);
%		str = 'coeff %-f %-f %-f %-f s2_coeff %-f %-f %-f %-f\n';
	otherwise
		error('get_model');
	end

	% amplitude (range) model
%	fD = @(c,k,D0,x,x0,Ur,r0,D1) D0.*exp(fr(c,k,x,x0,Ur,r0,D1).*(x-x0));

	% phase model
%	fk = fr;
	fz = @(c,k,D,phi0,x,x0,Ur,r0,D1) D.*exp(1i*(phi0-fk(c,k,x,x0,Ur,r0,D1)*(x-x0)));

	% set
	obj.c0 = c0;
%	obj.fr = fr;
	obj.fD = fD;
	obj.fk = @fk;
	obj.fz = fz;
end

% wave number (phase difference)
% power law in flow velocity
function fk = fk(c,k,x,x0,Ur,r0,D1)
	fk = -(c(1)./(x-x0) + c(2)*Ur.^c(3));
end

% amplitude models

% quadratic damping depending on river flow velocity
function [D r re] = fD1(c,k,D0,x,x0,Ur,r0,D1)
	c(3)  = 2;
	[D r re] = fD2(c,k,D0,x,x0,Ur,r0,D1)
end

% power law depending on flow velocity
function [D r re] = fD2(c,k,D0,x,x0,Ur,r0,D1)
	r = -( c(1)./(x-x0) ...
	       + c(2)*Ur.^c(3) ...
	     );
	D = D0.*exp(r*(x-x0));
	re = r;
end

% power law depending on flow velocity and self damping
% TODO, this can be done better, some theoretic support for coefficients
% TODO D1 must be D10
function [D r re] = fD3(c,k,D0,x,x0,Ur,r0,D1)
	if (1 == k) % major constituent with self damping
		r  = -(    c(1)./(x-x0) ...
                        + c(2)*Ur.^c(3) ...
			+ c(4)*r0.^2./sqrt(1+2*Ur) );
		D  = D0.*exp(r*(x-x0));
		re = r;
	else % secondary constituent including overtide generation
		r = -(    c(1)./(x-x0) ...
                        + c(2)*Ur.^c(3) );
		D = D0.*exp(r*(x-x0)) + c(4)*D1.^2.*Ur;
		% effective r
		re = 1./(x-x0)*log(D./D0);
	end
end

