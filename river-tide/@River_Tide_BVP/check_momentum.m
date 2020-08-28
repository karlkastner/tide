% Thu 10 Oct 11:36:29 PST 2019
% Karl Kastner, Berlin
%
%% compute residual for the momentum equation
%% dQ/dt + d/dx (Q^2/A) = - g A dz/dx - g A dw/dx - cd w Q|Q|/A^2
function [rmse,res] = check_momentum(obj)
	% TODO use flags for 1/h and advective acceleration
	for idx=1:obj.nc
		x   = obj.x(cdx);
		Q0  = obj.Q(0,cdx);
		zt  = obj.z(1,cdx);
		Qt  = obj.Q(1,cdx);
		w   = obj.width(cdx,x);
		h0  = obj.h0(cdx,x);
		cD  = obj.cd(cdx,x);

		[rmse(cdx),res{cdx}] = check_momentum@River_Tide(obj,x,w,zt,Qt,obj.omega);
	end
end % check_momentum

