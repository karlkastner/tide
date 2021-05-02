% Mon 20 Apr 15:35:09 +08 2020
% Karl Kastner, Berlin
%
% initialize the hydrodynamic solver
%
function y0 = initial_value(obj,x)
	x = cvec(x);

	% TODO, no magic numbers
	Q0_dummy = -10;
	h0_dummy =  10;
	z0_min   =   0;

	W  = obj.width(x);
	S  = obj.D1_dx(x)*obj.zb(x);
	Cz = drag2chezy(obj.cd(x, h0_dummy));
	h0 = normal_flow_depth(   Q0_dummy ...
				, W ...
				, Cz ...
				, S ...
			      );

	h0(~isfinite(h0)) = 0;
	h0 = max(h0,obj.rt.opt.hmin);

	z0 = obj.zb(x) + h0;
	z0 = max(z0,z0_min);

	%nf = sum(obj.opt.oflag);
	nf = sum(obj.rt.opt.oflag);
	e = 1;
	%if (~obj.opt.dischargeisvariable)
	if (~obj.rt.opt.dischargeisvariable)
		y0 = [z0; sqrt(eps)*randn(nf*length(x),1)];
	else
		%y0 = [z0; Q0_dummy; sqrt(eps)*randn(length(x)*(obj.neq-1),1)];
		%y0 = [z0+h0_dummy*randn(size(z0)); Q0_dummy*(1 + randn()); 1e-3*randn(length(x)*(obj.neq-1),1)];
		y0 = [z0; Q0_dummy; e*[randn(nf*length(x),1) + 1i*randn(nf*length(x),1)] ];
	end
end % initial_value

