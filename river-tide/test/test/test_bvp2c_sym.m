% Sun 29 Oct 20:47:55 CET 2017
function [x, y, dydx] = test_bvp2c
	%nx = 4;
	%L = nx-1;
	nx = 128;
	L  = 10;
	X = [0, L];
	opt = struct();
	opt.nx = nx;
	f0 = 2;
	ci = -1;
	r  =  0;
	rfun = @(x) r*x/L;
	

	neq=2;
	[x, y, cflag, dydx] = bvp2c(@odefun,@bcfun,X,opt);

	figure(1)
	clf();
	plot(x,[reshape(y,[],neq)]);%,(f0-1)*exp(-1*x)+rfun(x)]);

% simplify(dydx(1),'ignoreanalyticconstraints',true)

% -> expand already exp and eigenvalues as linear functions of x
	
	function [f p q] = bcfun(x,y0,cdx)
		switch (x)
		case {X(1)}
			f = f0;
			p = [1 0];
		case {X(end)}
			f = r;
			%p = [0 1];
			p = [1, 0];
		otherwise
			error('here')	
		end
		q = [1 1];
	end	

	function c = odefun(x,y)
		if (0 == nargin())
			c = zeros(0,0,neq);
			return;
		end
		for jdx=1:(neq)
		for idx=1:length(x)
			%rhs = 1;
			%r = 0;
			%r  = sym('r');
			%ci  = sym('-c0^2'); %+c1*x');
			%ci = sym(['c',num2str(idx)]);
			%eval(['syms c',num2str(idx)])
			%eval(['ci = c',num2str(idx)])
			c(idx,1:4,jdx)= [1, 0, ci, rfun(x(idx))];
		end
		end
	end
end

