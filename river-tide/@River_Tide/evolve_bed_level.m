% Fri  7 Aug 19:02:20 +08 2020
function [t,zb,zi] = evolve_bed_level(obj) %,zb0)
	% TODO, this is a workaround
	% evaluation of zb requires a mesh
	% but meshing requires a call to "solve"
	% meshing should be done in "init"

	obj.init();
	y0 = obj.solve();
	obj.opt.ifun = @(x) obj.out.ypm;

	y0 = obj.solve();

	zi     = obj.z_;
	x      = obj.x;
	zb0    = obj.zb(mid(x));	
	[t,zb] = obj.morsolver.solve(@obj.dzb_dt, zb0);
end

