% Thu 10 Oct 11:35:49 PST 2019
% Karl Kastner, Berlin

function [rmse, res] = check_continuity(obj,x,w,zt,Qt,omega)	
	% continuity : dA/dt + dQ/dx = 0
	% obj.D1_dx(cdx,x)*Q1;
	k     = (1:size(zt,2));
	res   = 1i*omega*w.*zt.*k + derivative1(x,Qt);
	rmse  = rms(res);
end % check_continuity

