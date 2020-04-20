% Mon 20 Apr 15:35:09 +08 2020
function y0 = ifun(obj,x)
	x = cvec(x);
	h0_dummy = 10;
	W = obj.width(x);
	S = obj.D1_dx(x)*obj.zb(x);
	cz = drag2chezy(obj.cd(x,h0_dummy));
	h0 =	normal_flow_depth(   obj.Q0_ ...
				   , W ...
				   , cz ...
				   , S ...
				  );

	h0(~isfinite(h0) | abs(S)<1e-6 | cz<1) = 0;
	% TODO this fails for right to left
	% h0 = max(obj.z0_downstream-obj.zb(obj.xi(1)), h0);
	h0 = max(obj.z0_downstream+h0_dummy, h0);
 
	z0 = obj.zb(x) + h0;
	y0 = [z0, rand(length(x),obj.neq-1)];
end

