% Mon 11 Dec 15:45:31 CET 2017

g     = Constant.gravity;
Q0    = 1;
w0    = 500;
h0    = 16;
omega = 2*pi/86400;
%az1  = 1e-4;
az1   = 0.1;
%qt0   = az1*sqrt(g*h0);
nx    = 1024;
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

r0 = damping_modulus(Q0,w0,h0,cd0,omega,az1)
R = [];
for idx=1:length(sigma)
	for jdx=1:2
	bc(1).q = [1 bcqr(jdx)];
	bc(2).p = [1 0];
	bc(2).rhs = 0;
	bc(2).q = [1 1];
	
	sw(2) = sigma(idx);
	
	hfun  = @(x) geometry_smooth(x, x0, h0_, sh)
	wfun  = @(x) geometry_smooth(x, x0, w0_, sw);
	cdfun = @(x) cd0 + (cde-cd0)*normcdf(x,x0(3),sc*lambda);
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

	[x y cflag dydx l c] = bvp2c(@out.odefun,@out.bcfun,Xi,opt);
	xc = 0.5*(x(2:end)+x(1:end-1));
	cQ  = reshape(c,2,[])';
	cQ  = cQ/abs(az1);
	cz  = -1./(1i*omega)*1./(wfun(xc)*[1 1]).*cdiff(cQ)./(cdiff(xc)*[1 1]);
	R(idx,1) = abs(cz(1,2));
	h0=out.h0;
	w=out.w(out.x);
	Rx(idx,1) = max(select(abs(1./w.*(cdiff(w)./cdiff(x))),1:nx/2));
	%Rx(idx,1) = max(select(abs(1./h0.*(cdiff(h0)./cdiff(x))),1:nx/2));
	
	%	q1_  = out.q1;
	%	z0_ = out.z0;
	%	k_ = [];
	figure(1);
	ns = length(sigma);
	subplot(2,ns,ns*(jdx-1)+idx)
	cla();
	plot(out.x,abs(z1));
	ylabel('z_1/|z_1(0)|')
	hold on
	plot(xc,abs(cz));
	%lot(out.x,exp(-r0*out.x));

	figure(2);
	subplot(2,ns,ns*(jdx-1)+idx);
	cla();
	plot(out.x,abs(Q1));
	hold on
	plot(xc,abs(cQ));
	ylabel('Q_1/|z_1(0)|')

	figure(3);
	subplot(2,ns,ns*(jdx-1)+idx);
	cla();
	plot(out.x,out.w(out.x));
	
	figure(4);
	subplot(2,ns,ns*(jdx-1)+idx);
	cla();
	plot(out.x,out.h0);
	ylabel('h0');

	figure(5);
	subplot(2,ns,ns*(jdx-1)+idx);
	cla();
	plot(out.x,out.cd(out.x));
	ylabel('cd')
	end % jdx
end % idx

figure(100);
% width convergence ~ 1/2 1/w dw/dx
% depth convergence ~ 1/4 1/h dw/dx
% we change width, so 
plot(2*Rx,R)
xlabel('R (reflection coefficient)');
ylabel('1/h dh/dx');
% approximate values for the Kapuas
h    = 4;
dhdx = 16/3e5;
vline(1/h*dhdx);
%sigma/lambda,R);

