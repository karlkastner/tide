% 2016-07-12 15:06:16.443700194 +0200
% Karl Kastner, Berlin
%
%% fit the tidal range
function [s2, R2, neff, obj] = fit_range(obj,time,R,R0,Ur,dUr_dt,x,x0)

	opt    = obj.opt;
	fD     = obj.fD;

	R0     = double(R0);
	par0   = double(obj.c0);
	R      = double(R);
	Ur     = double(Ur);
	dUr_dt = double(dUr_dt);

		fdx =   isfinite(dUr_dt) & isfinite(Ur) ...
                      & isfinite(R(:))  ...
		      & isfinite(R0(:));

		% fit
		[par resn res] = lsqnonlin(@(par) (R(fdx) - fD(par,1,R0(fdx),x,x0,Ur(fdx),R0(fdx),[])),par0,[],[],opt);

		% compute the hessian
%		H         = hessian(@(par) sum((D(fdx,idx)-fD(par,idx,D0(fdx,idx),x,x0,Ur(fdx),r0(fdx),D(fdx,1))).^2),par');
%		s2(idx)   = mean(res.^2);
%		R2(idx)   = 1 - s2(idx)/nanvar(D(fdx,idx));

		% effective sample size
%		rho  = nancorr(res(1:end-1),res(2:end));
%		rho  = rho(1,2);
%		f_   = (1+rho)/(1-rho);
%		neff(idx) = sum(fdx)/f_;
	%	s2_  = f_^2*s2/sum(fdx);
%		pars(:,idx) = s2(idx)/sqrt(neff(idx)-1)*diag(inv(H))';
%		par_(:,idx) = par;
%
%	end % for idx

	obj.c.range  = par;
%	obj.cs.amplitude = pars;
end

