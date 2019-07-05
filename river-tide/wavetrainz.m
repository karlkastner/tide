% Wed  1 Nov 11:40:08 CET 2017
% Karl Kastner, Berlin
%% determine river tide by iterated integration of the surface elevation
function [x, z1, q1, Q1, Z, QQ, F] = wavetrainz(cfun,X,z10,omega,hfun,wfun,opt)

	nout = nargout();

	if (nargin() < 7)
		opt = struct();
	end
	if (~isfield(opt,'nx'))
		nx = 1024;
	else
		nx = opt.nx;
	end
	if (~isfield(opt,'k'))
		k = 10;
	else
		k = opt.k;
	end
	
	% discretise domain
	dx = (X(2)-X(1))/(nx-1);
	x  = X(1) + dx*(0:nx-1)';
	w  = wfun(x);
	
	% initial estimate
	% TODO better estimate
	Q1 = ones(size(x));

	[Q1 cflag] = picard(@wavetrainz_,Q1);

	if (nout>4)
		Z = fliplr(cumsum(F.').');
		Z = scale*F;
		% TODO this misses the integration constant
		QQ = 1i*omega*cumintR(Z,dx);
	end

function Q1 = wavetrainz_(Q1)

	
	% coefficients of the wave equation
	[c ] = cfun(x);
	c = conj(c);
	
	% convergence term
	dcdx = cdiff(c)./dx;
	
	% reflection/transmission coefficient
	R   = +1/4*dcdx./c;
	
	% integral of the transmitted and damped travelling wave
	iL   = cumintL(-sqrt(-c) + R, dx);
	iL_  = cumintL(-sqrt(-c)    , dx);
	iR   = cumintR(-sqrt(-c) - R, dx);
	iR_  = cumintR(-sqrt(-c)    , dx);
	
	f0 = 1;
	
	% incoming wave
	f = f0*exp(iL);
	f_ = f0*((c(1)./c).^(-1/4)).*exp(iL_);

	if (nout>4)
		F  = zeros(nx,k);
		F_ = zeros(nx,k);
		F(:,1) = f;
		F_(:,1) = f_;
	end

	fs  = f;
	fs_ = f_;
	idx = 1;
	fs_old = 0;
	fsold = 0;
	while (idx<k)
		idx=idx+1;
		% waves reflected an odd number of times
		f = exp(iR).*cumintR( R.*f.*exp(-iR),dx);
		f_ = ((1./c).^(-1/4)).*exp(iR_) .* ...
                     cumintR(R.*f_.*(c).^(-1/4).*exp(-iR_),dx);
		fs = fs+f;
		fs_ = fs_+f_;
		%f_ = ((c(end)./c).^(1/4)).*exp(iR_) .* ...
                %     cumintR(R.*f_.*(c./c(end)).^(1/4).*exp(-iR_),dx);
	if (nout>4)
		F_(:,idx) = f_;
		F(:,idx)  = f;
	end


		if (idx>k)
			break;
		end
		idx = idx+1;
	
		% waves reflected an even number of times
		f = exp(iL).*cumintL(-R.*f.*exp(-iL),dx);
		f_ = (1./c).^(-1/4).*exp(iL_) .* ...
		     cumintL(-R.*f_.*(c).^(-1/4).*exp(-iL_),dx);
		%f_ = (c(1)./c).^(1/4).*exp(iL_) .* ...
		%     cumintL(-R.*f_.*(c./c(1)).^(1/4).*exp(-iL_),dx);
		fs = fs+f;
		fs_ = fs_+f_;
		if (nout>4)
			F(:,idx) = f;
			F_(:,idx) = f_;
		end % if nargout

		if (max(abs(fs - fsold)) < sqrt(eps))
			idx
			break;
		end
		fsold = fs;
		fs_old = fs_;
	end % while
	% norm(F-F_)
	
	z1 = fs;

	% apply initial value
	scale = z10/z1(1);
	z1 = scale*z1;
	
	% compute discharge
%	q = 1i*omega*cumintR(z,dx);
	q1 = z2q(dx,z1,hfun,omega);
	Q1 = w.*q1;

end % wavetrainz_

end % wavetrainz

