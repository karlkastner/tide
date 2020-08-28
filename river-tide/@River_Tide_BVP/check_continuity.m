% Thu 10 Oct 11:35:49 PST 2019
function [rmse,res] = check_continuity(obj)
	for cdx=1:length(obj.nx)
		x     = obj.x(cdx);
		omega = obj.omega;
		w     = obj.width(cdx,x);
		Q1    = obj.Q(cdx,1);
		z1    = obj.z(cdx,1);

		[rmse(cdx),res{cdx}] = check_continuity@River_Tide(obj,x,w,zt,Qt,omega);
	end
end % check_continuity

