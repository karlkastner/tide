% Wed 11 Oct 10:18:54 CEST 2017
% Karl Kastner, Berlin
%
%% postprocess hydrodynamic solver
%
function postprocess(obj,cdx,x,y,yc)	
	% extract unknowns
	nx         = length(x);
	nxc	   = nx-1;
	xc         = mid(x);
	w          = obj.width(cdx,x);
	wc         = mid(w);

	[z0,   Q0,  Qt] = obj.extract(x,y);
	[z0c, Q0c, Qtc] = obj.extract(xc,yc);

	% stack
	obj.out(cdx).Q  = [Q0*ones(size(x)), Qt];
	obj.out(cdx).z  = [z0, obj.discharge2level(x,Qt,w)];

	obj.out(cdx).Qc = [Q0*ones(nxc,1), Qtc];
	obj.out(cdx).zc = [z0c, obj.discharge2level(xc,Qtc,wc)];
	
	% check result
	cflag = 0;
	if (any(obj.h0(cdx,x)<=0))
		cflag = -2;
		warning('negative water depth');
	end

	if (any(~isfinite(obj.out(cdx).z)))
		cflag = -1;
		warning('solution is not finite');
	end

	obj.out.cflag_ = cflag;
end

