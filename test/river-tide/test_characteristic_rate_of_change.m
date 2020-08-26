% Tue 12 Dec 14:19:11 CET 2017

function test
g     = Constant.gravity;
Q0    = 1;
w0    = 500;
h0    = 16;
omega = 2*pi/86400;
%az1  = 1e-4;
az1   = 0.1;
%qt0   = az1*sqrt(g*h0);
nx    = 4*1024;
lambda = 2*pi*sqrt(g*h0)/omega;

%Xi = [0, L];
%L     = 20e6;
L =   20*lambda;
Xi    = [0,L];

sigma = [0,1];

	
	%h0_   = h0*[1 (1/16) 0.001];
	h0_   = h0*[1 1 0.001];
	%h0_   = h0*[1 1/16 1/16];
	sh    = [0 0.5 1]*lambda;
	sc    = 1;
	x0    = Xi(1)+[[0,1/4]*lambda,L-2*sc*lambda];


	% h01/p^4;
	w0_   = w0*[1, 1/4, 1];
	sw    = lambda*[0,1,1];
	sigma = lambda*[1/2,1/4,1/8,1/16];
	%lambda*[1e-2 1e-1 1];
%	ch   = -4*log(p)/(0.5*L);

	cd0   = 0.0;
	cde   = 0.1;

opt.nx = 1024;

bc  = struct();
bcqr = [1 0];
iflag = false;

flipflag = true;
%flipflag = false;

r0 = damping_modulus(Q0,w0,h0,cd0,omega,az1)
R = [];

%	for jdx=1:2
	jdx=1;
	bc(1).q = [1 bcqr(jdx)];
	bc(1).p = [1 0];

	bc(2).p = [1 0];
	bc(2).rhs = 0;
	bc(2).q = [1 1];
	
%	sw(2) = sigma(1);
	
%	hfun  = @(x) geometry_smooth(x, x0, h0_, sh)
%	wfun  = @(x) geometry_smooth(x, x0, w0_, sw);
	hfun  = @(x) h0*ones(size(x));

	w0    = 1e5;
	we    = 500;
%	w0 = we;

	
	s = 0.5;
	%s = 12;
	s = 10;
	p = 1;
	if (flipflag)
	%	wfun  = @(x) geometry_smooth(Xi(2)-x, [0 0 4 (Xi(2)/lambda-8)]*lambda, [w0 we w0 we], [0 2*lambda s*lambda 4*lambda])
		wfun  = @(x) we + (w0-we)*exp(-(s/lambda*(Xi(2)-x)).^p);
	else
		wfun  = @(x) we + (w0-we)*exp(-(s/lambda*x).^p);
	%	wfun  = @(x) geometry_smooth(x, [0 0 4 (Xi(2)/lambda-8)]*lambda, [w0 we w0 we], [0 2*lambda s*lambda 4*lambda])
	end
	%wfun  = @(x) geometry_smooth(x, [0 0 4 (Xi(2)/lambda-8)]*lambda, [w0 we w0 we], [0 2*lambda 0.05*lambda 4*lambda])
	%wfun  = @(x) geometry_smooth(x, [0 4]*lambda, [w0 we], [0 0.05*lambda])
	cdfun = @(x) 2e-3*ones(size(x));
	%cdfun = @(x) cd0 + (cde-cd0)*normcdf(x,x0(3),sc*lambda);
	%Xi(2)-sc*lambda,sc*lambda);

%	if (~iflag)
	out = RT(  'fun.zb', @(x) -hfun(x) ...
			, 'cd', cdfun ...
			, 'w', wfun ...
			, 'omega', omega ...
			, 'bc', bc ...
			, 'opt', opt ...
			... %, 'solver', @bvp2fdm ...
			, 'solver', @bvp2c ...
			... %, 'flag', flag ...
	);
		out.init([0,az1], Q0, Xi);		% TODO q0*w
		iflag = true;
%	else
%		out.bc = bc;
%		out.Q1 = [];
%	end
	out.solve();
	%	x_  = out.x;
	z1 = out.z1/abs(az1);
	Q1 = out.Q1/abs(az1);

	[x y cflag dydx0 l c] = bvp2c(@out.odefun,@out.bcfun,Xi,opt);
	xc  = 0.5*(x(2:end)+x(1:end-1));
	yc  = reshape(c,2,[])'; %/az1;
	%yc = conj(yc);

	% I completely do not know why the values have to be conugated here
%	yc = conj(yc);
	l = conj(l);
%	l = [l(:,2),l(:,1)];

%colnorm([real(y) real(yc)])
%colnorm([i(yc) imag(yc)])
%pause

%	l = [l(:,2),l(:,1)];
	%dydx  = diff(yc)./(diff(xc)*[1 1]);
	%dydx(end+1,:) = NaN;
	dydx  = cdiff(yc)./(cdiff(xc)*[1 1]);
	dydx_ = ode2characteristic(xc,l,yc);

	figure(1)
	clf();
	subplot(2,2,1)
	plot(x,abs(z1))
	subplot(2,2,2)
	plot(x,abs(Q1))
	subplot(2,2,3)
	plot(x,real(y))
	hold on
	plot(xc,real(yc))
	subplot(2,2,4)
	plot(x,imag(y))
	hold on
	plot(xc,imag(yc))

%	subplot(2,2,4)
%	plot(x,wfun(x))
	%plot(xc,abs(y))
%$	subpl

	n=4;
	if (flipflag)
		fdx = xc>Xi(2)-n*lambda;
	else
		fdx = xc<n*lambda;
	end
	res = dydx(fdx,:)-dydx_(fdx,:);
	[colnorm([real(res(:,1)) (imag(res(:,1)))]) colnorm([real(res(:,2)) (imag(res(:,2)))])]

	n = 5;
	if (flipflag)
		xlim_ = [Xi(2)-n*lambda,Xi(2)];
	else
		xlim_ = [Xi(1),Xi(1)+n*lambda];
	end
	

	% make rates relative
%	r = r./abs(y);
%	dydx = dydx./abs(y);

	figure(2);
	clf();
	for idx=1:2
		subplot(2,2,idx);
		plot(xc,real([dydx(:,idx) dydx_(:,idx)]));
		hold on
		plot(xc,real(l(:,idx).*yc(:,idx)),'--');
		xlim([xlim_]);
		%hold on
		subplot(2,2,2+idx);
		plot(xc,imag([dydx(:,idx) dydx_(:,idx)]));
		hold on
		plot(xc,imag(l(:,idx).*yc(:,idx)),'--');
		%plot(xc,imag(l.*yc),'--');
		xlim([xlim_]);
	end

	dldx = cdiff(l)./(cdiff(xc)*[1 1]);
	figure(3);
	clf();
	for idx=1:2
		subplot(2,2,idx);
		d = l(:,2)-l(:,1);
		%plot(xc,[real(l(:,idx)) imag(l(:,idx))]);
		plot(xc,real([l(:,idx)  (1./d.*dldx(:,idx))]))
		xlim([xlim_]);
		%hold on
		subplot(2,2,2+idx);
		%plot(xc,[real(dldx(:,idx)) imag(dldx(:,idx))]);
		%plot(xc,[real(1./d.*dldx(:,idx)) imag(1./d.*dldx(:,idx))]);
		plot(xc,imag([l(:,idx)  (1./d.*dldx(:,idx))]))
		xlim([xlim_]);
	end

	figure(4);
	clf();
	plot(x,wfun(x))

	figure(6);
	if (flipflag)
		yc_ = bsxfun(@plus, yc(end,:), [cumintR(dydx_(:,1),xc),cumintR(dydx_(:,2),xc)]);
	else
		yc_ = bsxfun(@plus, yc(1,:), [cumintL(dydx_(:,1),xc),cumintL(dydx_(:,2),xc)]);
	end
%	subplot(2,2,1);
%	plot(xc,[real(yc_(:,1)));
%	subplot(2,2,2);
%	plot(xc,yc_(:,2));

	for idx=1:2
		subplot(2,2,idx);
		plot(xc,real([yc(:,idx)  yc_(:,idx)]))
		xlim([xlim_]);

		subplot(2,2,2+idx);
		plot(xc,imag([yc(:,idx)  yc_(:,idx)]))
	end
end


