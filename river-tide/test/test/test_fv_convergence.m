% 2017-11-14 23:45:17.885372831 +0100

	% convergence test
	N = 8*(2.^(1:7))+1;
	R = [];
	for idx=1:nlim
	figure(100+idx);
	clf
	for jdx=length(N):-1:1
		nx = N(jdx);
		switch (limiter_C{idx})
		case {'lax_friedrich'}
			fv = Lax_Friedrich();
		otherwise
			fv         = Reconstruct_Average_Evolve();
			%FL = Flux_Limiter();
			fv.limiter = @(varargin) Flux_Limiter.(limiter_C{idx})(varargin{:});
		end
		fv.bcfun{1} = @SWE.bc_nonreflecting;
		fv.bcfun{2} = @SWE.bc_nonreflecting;
		fv.icfun    = ic;
		fv.pde      = SWE();
		fv.init([-L/2,L/2],nx);
		[T Y] = fv.solve(Ti);
		h = Y(1:nx,end);
		figure(100+idx);	
		plot(fv.x,h);
		hold on
		if (length(N)==jdx)
			xref       = fv.x;
			href       = h;
			R(jdx,idx) = NaN;
		else
			h = interp1(fv.x(1+k:end-k),h(1+k:end-k),xref,'spline');
			R(jdx,idx) = norm(h-href);
		end
	end % for jdx
	end % for idx
	figure(1000)
	loglog(N,R);
	legend(limiter_C{1:nlim})


