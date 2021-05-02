% Mon  5 Oct 13:31:46 +08 2020
function x0 = mg_step(obj,x0)
	% the up-of the mg-v-cycle is not part here
	% as local bed level changes near channel ends might be too large as long
	% as the low frequency components have not yet been computed
	%
%	[A,b] = fun(x0);
%	x0    = A \ b;
	if (any(obj.nxc>3))
		x1    = obj.mg_restrict(x0);
		x1_   = obj.lower.mg_step(x1);
		x0_   = x0 + obj.mg_interpolate(x1_-x1);
	end % if
% TODO this has to be the single bed-level step
	[A,b] = fun(x0,x0_);
	x0    = A \ b;
end % function mg_step

