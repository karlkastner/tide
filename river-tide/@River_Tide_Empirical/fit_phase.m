% 2016-07-12 15:06:16.443700194 +0200
% Karl Kastner, Berlin
%
%% fit the phase of the oscillatory components
function [s2, R2, neff, obj] = fit_phase(obj,time,z,z0,Ur,dUr_dt,x,x0,r0)

	D0 = abs(z0);
	D  = abs(z);	
	phi0 = angle(z0);
	phi  = angle(z);

	delta_phi = phi - phi0;
	delta_phi = phasewrap(delta_phi);

	opt = obj.opt;
	fz  = obj.fz;

%	[fD fr fz fk void0 par0] = rt_model(model);

	% cast, as lsqnonlin requires double input
	D0     = double(D0);
	par0   = double(obj.c0);
	D      = double(D);
	Ur     = double(Ur);
	r0     = double(r0);
	dUr_dt = double(dUr_dt);

	% for each constituent
	for idx=1:size(D0,2)
		fdx = isfinite(dUr_dt) & isfinite(Ur) ...
		      & isfinite(D(:,idx)) & isfinite(D0(:,idx)) ...
		      & isfinite(D(:,1));
		% fit
		[par resn res]  = lsqnonlin(@(par) objective(z(fdx,idx) - fz(par,idx,D(fdx,idx),phi0(fdx,idx),x,x0,Ur(fdx),r0(fdx),D(fdx,1))),par0,[],[],opt);
		% compute the hessian
		H    = hessian(@(par) sum((z(fdx,idx)-fz(par,idx,D(fdx,idx),phi0(fdx,idx),x,x0,Ur(fdx),r0(fdx),D(fdx,1))).^2),par');
		s2(idx)   = rms(res).^2;
		R2(idx)   = 1 - s2(idx)/nanvar(z(fdx,idx));
		% effective sample size
		rho  = nancorr(res(1:end-1),res(2:end));
		rho  = rho(1,2);
		f_   = (1+rho)/(1-rho);
		neff(idx) = sum(fdx)/f_;
	%	s2_  = f_^2*s2/sum(fdx);
		pars(:,idx) = s2(idx)/sqrt(neff(idx)-1)*diag(inv(H))';
		par_(:,idx) = par;
	end % for idx
	%a    = nanautocorr(res,200);
	%rho_ = ar1delay(a,1)
	%f__ = sqrt((1+rho_)/(1-rho_));
	%[sum(fdx),f_.^2,f__.^2]
	obj.c.phase  = par_;
	obj.cs.phase = pars;
end

function res = objective(res)
	res = [real(res);imag(res)];
end

