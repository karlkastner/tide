% Fri  7 Aug 19:02:20 +08 2020
% TODO use the non-linear multigrid method to evolve over larger time spans
function [t, zb, zi] = evolve_bed_level(obj)
	% TODO this is a workaround
	% evaluation of zb requires a mesh
	% but meshing requires a call to "solve"
	% meshing should be done in "init"

	obj.init();
	obj.solve();

	zi = [];

	% store solution for reuse as initial condition
	for cdx=1:obj.nc
		obj.tmp(cdx).ypm = repmat(obj.odesolver.out(cdx).ypm,1,obj.season.iorder);
		if (nargout()>2)
			zi_ = obj.out(cdx).z;
			zi(1:size(zi_,1),end+1:end+size(zi_,2)) = zi_;
		end
	end

	zb0 = [];
	ni = 0;
	for cdx=1:obj.nc
		% TODO fetch xc
		x      = obj.odesolver.out(cdx).x;
		nxc    = length(x)-1;
		zb0(ni+(1:nxc),1) = obj.zb(cdx,mid(x));
		ni = ni+nxc;
	end

	[t,zb] = obj.morsolver.solve(@obj.dzb_dt, zb0);

	obj.evolution.t  = t;
	obj.evolution.zb = zb;
end % evolve bed level

