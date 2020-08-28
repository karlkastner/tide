% Mon 20 Apr 15:35:09 +08 2020
function y0 = ifun(obj,cdx,x)
	x = cvec(x);
	% TODO, no magic numbers
	Q0_dummy = -10;
	h0_dummy =  10;
	z0_dummy =   0;

	W  = obj.width(cdx,x);
	S  = obj.D1_dx(cdx,x)*obj.zb(cdx,x);
	cz = drag2chezy(obj.cd(cdx,x,h0_dummy));
	h0 =	normal_flow_depth(   Q0_dummy ...%obj.Q0_ ...
				   , W ...
				   , cz ...
				   , S ...
				  );

	h0(~isfinite(h0) | abs(S)<1e-6 | cz<1) = 0;
	% TODO this fails for right to left
	% h0 = max(obj.z0_downstream-obj.zb(obj.xi(1)), h0);
	%h0 = max(obj.z0_downstream+h0_dummy, h0);
	h0 = max(z0_dummy+h0_dummy, h0);
 
	z0 = obj.zb(cdx,x) + h0;

	if (~obj.opt.dischargeisvariable)
		y0 = [z0; rand(length(x),obj.neq-1)];
	else
		%y0 = [z0; Q0_dummy; sqrt(eps)*randn(length(x)*(obj.neq-1),1)];
		y0 = [z0+h0_dummy*randn(size(z0)); Q0_dummy*(1 + randn()); 1e-3*randn(length(x)*(obj.neq-1),1)];
	end
end

