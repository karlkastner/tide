% Mon  5 Oct 13:32:44 +08 2020
function mg_prepare(obj,x,nx)
	if (nargin()>1)
		obj.x  = x;
		obj.nx = nx;
	else
		x = obj.x;
	end
	if (any(obj.nxc>3))
		obj.lower = River_Tide_BVP(obj);
		% for each channel
		for idx=1:obj.nc
			if (obj.nxc(cdx)>3)
				% note, this is also different to the standard multigrid
				% as the bed level is sampled at the segment centres
				%xp(TODO)  = x(TODO) 2:2:end-1);
				% TODO this only holds for nxc
				nxp(cdx)  = 0.5*(nx+1)-1;
			else
				xp(TODO) = x(TODO)
				nxp      = nx(cdx);
			end
		end
		% recurse
		% TODO compute indices
		obj.lower.mg_prepare(xp,nxp);
	end % if any
end % mg_prepare


