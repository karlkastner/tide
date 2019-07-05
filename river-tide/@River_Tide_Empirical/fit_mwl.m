% Mon 11 Jul 15:42:54 CEST 2016
% Karl Kastner, Berlin
%
%% fit the tidally averaged water level
% h : low pass filtered depth (can have offset)
function [rms_, R2, obj] = fit_mean_level(obj,h,R0,Qr);
	Qr    = double(Qr);
	R0    = double(R0);

	Qr   = abs(Qr);
	fdx  = isfinite(h) & isfinite(R0) & isfinite(Qr);
	n    = sum(fdx);

	% set up regression matrix
	X    = [ones(n,1),Qr(fdx).^(2/3), abs(R0(fdx)).^2./(Qr(fdx).^(4/3))];

	% linear regression
	[Q, R] = qr(X,0);
	c      = R \ (Q'*h(fdx));

	opt = obj.opt;

	% non linear regression
	if (obj.nlflag)
		[c, resn, res] = lsqnonlin(@(c) h(fdx) - (c(1) + c(2)*(Qr(fdx).^2 + c(3)*abs(R0(fdx)).^2).^(1/3)),c,[],[],opt);
	else
		res = X*coeff - h(fdx);
	end
	
	rms_ = rms(res);
	R2   = 1 - rms_^2/var(h(fdx));

	obj.c.mwl = c;
end

