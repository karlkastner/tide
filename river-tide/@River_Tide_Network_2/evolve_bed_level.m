% Fri  7 Aug 19:02:20 +08 2020
% TODO use the non-linear multigrid method to evolve over larger time spans
function [t,zb,zi] = evolve_bed_level(obj)
	% TODO this is a workaround
	% evaluation of zb requires a mesh
	% but meshing requires a call to "solve"
	% meshing should be done in "init"
	% TODO bvp2c should be an object, storing the mesh

	obj.init();
	y0 = obj.solve();

	zi = [];

	% store solution for reuse as initial condition
	for idx=1:length(obj.rt)
		%obj.rt(idx).opt.ifun = @(x) obj.rt(idx).out.ypm;
		obj.tmp(idx).ypm = repmat(obj.rt(idx).out.ypm,1,obj.season.iorder);
		if (nargout()>2)
			zi_ = obj.rt(idx).z_;
			zi(1:size(zi_,1),end+1:end+size(zi_,2)) = zi_;
		end
	end

%	obj.opt.ifun = @(x) y0;
%	zi     = obj.z_;
%	nci = obj.nci;
	zb0 = [];%zeros(nci(end)-1,1);
	ni = 0;
	for idx=1:length(obj.rt)
		% TODO fetch xc
		x      = obj.rt(idx).x;
		nxc    = length(x)-1;
		%zb0(nci(idx):nci(idx+1)-1) = obj.rt(idx).zb(mid(x));	
		zb0(ni+(1:nxc),1) = obj.rt(idx).zb(mid(x));
		ni = ni+nxc;
	end
	[t,zb] = obj.rt(1).morsolver.solve(@obj.dzb_dt, zb0);

	obj.evolution.t  = t;
	obj.evolution.zb = zb;
end

