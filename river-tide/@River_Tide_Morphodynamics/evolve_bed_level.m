% Fri  7 Aug 19:02:20 +08 2020
function [t,zb,zi] = evolve_bed_level(obj,zb0,T,nt)
	obj.init();
	obj.solve();
	zi = obj.z_;
%	zb         = obj.zb(mid(obj.x));
%	obj.fun.zb = @(x) interp1(mid(obj.x),zb,x,'linear');
	obj.init();
	obj.solve();
	x   = obj.x;
	zb0 = obj.zb(mid(x));	
	norder = obj.norder;
	[t,zb] = ode_adams_bashforth(@(t,zb) obj.dzb_dt(t,zb),zb0,T,nt,norder,obj.ks);
end

