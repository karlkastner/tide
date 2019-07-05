% 2016-07-12 15:06:16.443700194 +0200
% Karl Kastner, Berlin
%
%% fit the oscillatory components
function [s2, R2, neff, obj] = fit_amplitude(obj,time,z,z0,Ur,dUr_dt,x,x0,r0)
	D0     = abs(z0);
	D      = abs(z);

	opt    = obj.opt;
	fD     = obj.fD;

	D0     = double(D0);
	par0   = double(obj.c0);
	D      = double(D);
	Ur     = double(Ur);
	r0     = double(r0);
	dUr_dt = double(dUr_dt);
	% for each constituent (or species)
	for idx=1:size(D0,2)
		fdx =   isfinite(dUr_dt) & isfinite(Ur) ...
                      & isfinite(D(:,idx))  ...
		      & isfinite(D0(:,idx)) ...
		      & isfinite(r0) ...
		      & isfinite(D(:,1));

		% fit
		[par resn res] = lsqnonlin(@(par) (D(fdx,idx) - fD(par,idx,D0(fdx,idx),x,x0,Ur(fdx),r0(fdx),D(fdx,1))),par0,[],[],opt);

		% compute the hessian
		H         = hessian(@(par) sum((D(fdx,idx)-fD(par,idx,D0(fdx,idx),x,x0,Ur(fdx),r0(fdx),D(fdx,1))).^2),par');
		s2(idx)   = mean(res.^2);
		R2(idx)   = 1 - s2(idx)/nanvar(D(fdx,idx));

		% effective sample size
		rho  = nancorr(res(1:end-1),res(2:end));
		rho  = rho(1,2);
		f_   = (1+rho)/(1-rho);
		neff(idx) = sum(fdx)/f_;
	%	s2_  = f_^2*s2/sum(fdx);
		pars(:,idx) = s2(idx)/sqrt(neff(idx)-1)*diag(inv(H))';
		par_(:,idx) = par;

	end % for idx

	obj.c.amplitude  = par_;
	obj.cs.amplitude = pars;
end

