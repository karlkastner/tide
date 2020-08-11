% Wed 11 Oct 10:18:54 CEST 2017
function postprocess(obj,x,y)	
	% extract unknowns
	nx    = length(x);
	obj.x = x;
	obj.opt.nx = length(x);
	y_    = reshape(y,nx,[]);

	switch (obj.opt.hmode)
	case {'matrix'}
		z0 = y_(:,1);
		Qt = y_(:,2:end);
	otherwise
		z0 = obj.tmp.z0;
		Qt = y_;
	end % switch

	% stack
	obj.Q_ = [obj.Q0_*ones(size(x)), Qt];
	obj.z_ = [z0, obj.discharge2level(Qt)];
	
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

