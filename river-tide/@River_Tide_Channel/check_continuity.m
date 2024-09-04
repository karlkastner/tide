% Thu 10 Oct 11:35:49 PST 2019
% Karl Kastner, Berlin
%
%% compute residual for the continuity equation
%% dA/dt + dQ/dx = Q_in
%%
function [rmse, res] = check_continuity(obj)
	x     = obj.x;
	w     = obj.width();
	Q1    = obj.discharge(1);
	z1    = obj.waterlevel(1);

	[rmse, res] = obj.rt.check_continuity(x,w,z1,Q1);
end

