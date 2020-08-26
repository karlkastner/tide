% Wed 11 Oct 10:18:54 CEST 2017
function postprocess(obj,x,y,yc)	
	% extract unknowns
	nx         = length(x);
	nxc	   = nx-1;
	obj.x      = x;
	xc         = mid(x);
	obj.opt.nx = length(x);

	[z0,   Q0, Qt] = obj.extract(x,y);
	[z0c, Q0c, Qtc] = obj.extract(xc,yc);

	
	obj.Q0_ = Q0;

	% stack
	obj.Q_  = [Q0*ones(size(x)), Qt];
	obj.z_  = [z0, obj.discharge2level(Qt)];

	obj.Qc_ = [Q0*ones(nxc,1), Qtc];
	obj.zc_ = [z0c, obj.discharge2level(Qtc,mid(x))];
	
	% check result
	cflag = 0;
	if (any(obj.h0(x)<=0))
		cflag = -2;
		warning('negative water depth');
	end
	if (any(~isfinite(obj.z_)))
		cflag = -1;
		warning('solution is not finite');
	end

	obj.out.cflag_ = cflag;
end

