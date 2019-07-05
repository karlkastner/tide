% Mon  4 Dec 14:10:05 CET 2017
% Karl Kastner, Berlin
%% two pass solution for the linearised wave equation, for surface elevation
function [x, z, q, k] = wave_twopassz(K2fun,Xi,omega,opt,hfun)
	flag = 0;
	% domain length
	L = Xi(2)-Xi(1);

	if (nargin()<4)
		opt = struct();
	end
	if (~isfield('opt','nx'))
		opt.nx = 1024;
	end

	opt.MaxStep = L/opt.nx;

	% initial value for k
	k_L = sqrt(-K2fun(L));

	% from L to zero: solve for k
	[xk k] = ode23(@kdot,[Xi(2),Xi(1)],k_L,opt);
	%[xk k] = ode23(@kdot,Xi,k_L,opt);
	% flip x back
	if (flag)
		xk = flipud(L-xk);
		k  = flipud(k);
	end

	dkdx  = cdiff(k)./cdiff(xk);
	kfun  = @(xi) interp1(xk,k,xi);
	dkfun = @(xi) interp1(xk,dkdx,xi);

	% from 0 to L: solve for q
	z0    = 1;
	[x, z] = ode23(@zdot,Xi,z0,opt);
	% resample k to x
	k     = kfun(x);	

	% determine discharge from surface elevation
	q = z2q(x,z,hfun,omega);

	% apply bc
	scale = 1/z(1);
	z = scale*z;
	q = scale*q;
	
	function kdot = kdot(x,k)
		% solve backward in space
		if (flag)
			x   = L-x;
		end
		% TODO invert sign if flag is set
		kdot = -(k.^2 + K2fun(x));
	end
	
	function zdot = zdot(x,z)
		k_    = kfun(x);
		zdot  = (k_+1./k_.*dkfun(x)).*z;
	end	
end % wave_twostep

