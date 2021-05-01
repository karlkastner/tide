% Fri  7 Aug 19:02:20 +08 2020
% Karl Kastner, Berlin
%
%% evolve the bed level of the tidal river network over time
%
% TODO use the non-linear multigrid method to evolve over larger time spans
function [t, zb] = evolve_bed_level(obj)

	obj.init();

	zb0 = [];
	ni  = 0;
	for cdx=1:obj.nc
		% TODO fetch xc
		x      = obj.hydrosolver.out(cdx).x;
		nxc    = length(x)-1;
		zb0(ni+(1:nxc),1) = obj.zb(cdx,mid(x));
		ni = ni+nxc;
	end

	[t,zb] = obj.morsolver.solve(@obj.dzb_dt, zb0);

	obj.evolution.t  = t;
	obj.evolution.zb = zb;
end % function evolve_bed_level

