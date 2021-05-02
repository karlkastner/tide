% Wed 11 Oct 10:18:54 CEST 2017
% Karl Kastner, Berlin
%
%% postprocess hydrodynamic solver output
%
function cflag = postprocess(obj,x,y,yc)	
	% obj.rt.hydrosolver.out(obj.id).x;
	obj.x = x;

	% extract unknowns
	nx         = length(x);
	nxc	   = nx-1;
	xc         = mid(x);
	w          = obj.width(x);
	wc         = mid(w);

	[z0,   Q0,  Qt] = obj.extract(x,y);
	[z0c, Q0c, Qtc] = obj.extract(xc,yc);

	% stack
	obj.Q  = [Q0*ones(size(x)), Qt];
	obj.z  = [z0, obj.rt.discharge2level(x,Qt,w)];

	obj.Qc = [Q0*ones(nxc,1), Qtc];
	obj.zc = [z0c, obj.rt.discharge2level(xc,Qtc,wc)];
	
	% check result
	cflag = 0;
	if (any(obj.h0(x)<=obj.rt.opt.hmin))
		cflag = -2;
		warning('negative water depth (zs<hmin)');
	end
%	if (min(h0)<=0)
%		warning('negative water depth')
%	end

	if (any(~isfinite(obj.z)))
		cflag = -1;
		warning('solution is not finite');
	end

	% TODO parent
	%obj.cflag = cflag;
end

