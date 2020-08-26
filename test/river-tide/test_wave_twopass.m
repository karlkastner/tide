% Tue  5 Dec 01:21:13 CET 2017
function out = test_wave_twopass()

g     = Constant.gravity;
Q0    = 0;
Q0    = 1;
w0    = 500;
q0    = Q0/w0;
h0    = 16;
omega = 2*pi/86400;
L     = 20e6;
%az1  = 1e-4;
az1   = 1e-4;
qt0   = az1*sqrt(g*h0);
nx    = 1024;

Xi = [0, L];

% TODO rigid lid/constant depth mode for RT wave-equation solve

% rt:wrong
%leg = {'w2z'};
%leg = {'w2q'};
%leg  = {'rt'};
%leg  = {'wtz'};
%leg  = {'wtq'};
%leg = {'w2z','w2q','wtz','wtq'};
%leg = {'rt', 'wtq'} %,'wtq','w2q'};
%leg = {'rt', 'w2q'}
leg = {'rt','wtq','w2q'};
%leg = {'rt','w2q','w2z','wtz','wtq'};
%leg = {'rt','w2q','wtq'};
%leg = {'rt','w2q','wtq'};

lambda = 2*pi*sqrt(g*h0)/omega;

for idx=2
%1:1
if (1 == idx)
	% case 1: damping in uniform channel
	% constant depth
	hfun = @(x) h0*ones(size(x));
	cd    = 2e-3;
else
	ph1  = 0.5; % 2
	sh1  = 0.1;
	ph2  = 2; % 0.5
	sh2  = 0.1;	
	h01  = h0*ph1.^-4;
	h02  = h01*ph2.^-4;
	h0e  = min(h02,1);
	she   = 2;


	% h01/p^4;
	pw1   = 1;
	pw2   = 1;
	ww1   = 1;
	ww2   = 1;
	w01  = w0*pw1.^-2;
	w02  = w01*pw2.^-2;
%	ch   = -4*log(p)/(0.5*L);
	hfun = @(x) geometry_smooth(x,[0,1/4,2/4,3/4]*L, [h0, h01, h02, h0e],[0,sh1,sh2,she]*lambda);
	wfun = @(x) geometry_smooth(x,[0,1/4,2/4],[w0,w01,w02],[0,ww1,ww2]*lambda);

end
%	cdfun = @(x) cd.*(x<=0.5*L) + 1*(x>0.5*L)
%	cdfun    = @(x) fun(x,cd,1e-3);
	cd0 = 1e-5;
	cd0 = 1e-3;
	cd0 = 1e-1;
	%cd0 = 1;
	ce  = 0.2;
	cdfun = @(x) cd0 + ce*normcdf(x,3/4*L,she*lambda);
%	cdfun = @(x) 1e-3*normcdf(x,3/L,L/8);
%	cdfun = @(x) cd*ones(size(x));

	
	bc(2).rhs = 0;
	% TODO, this is a non-linear bc actually
	%bc(2).p = [(- omega.^2./(g.*hfun(L)) + 8/(3*pi).*1i.*omega.*cdfun(L)./g.*qt0./hfun(L).^3),0,-1];
	bc(2).p = [0 1 0];
	%bc(2).p = [0 1 0];
	% TODO flags

    opt=struct();
    %nx = 1024;
    %nk = 100;
    opt.obc = 1;
    opt.nx = nx;;
    opt.kmax  = nx;
    z0 = NaN(opt.nx,1);
    x = linspace(Xi(1),Xi(2),opt.nx)';
    z1 = [];
    k = [];
    q1 = [];
    k__ = NaN(nx,1);
    for jdx=1:length(leg)
	switch (leg{jdx})
	case {'rt'}
		out = RT(  'fun.zb', @(x) -hfun(x) ...
				, 'cd', cdfun ...
				, 'w', wfun ...
				, 'omega', omega ...
				, 'bc', bc ...
				, 'opt', opt ...
				... %, 'flag', flag ...
			);
		out.init([0,az1], Q0, Xi);		% TODO q0*w
		out.solve();
		x_  = out.x;
		z1_ = out.z1;
		q1_  = out.q1;
		z0_ = out.z0;
		k_ = [];
	case {'w2z'}
		[x_, z1_, q1_, k_] = wave_twopassz(@cfun,Xi,omega,opt,hfun);
	case {'w2q'}
		% [x_, z1_, q1_, k_] = wave_twopassq(@cfun,Xi,omega,opt);
		out2 = out();
		% reset
		out2.Q1 = [];
		out2.solver = @bvp2wavetwopass;
		out2.solve();
		x_  = out2.x;
		z1_ = out2.z1;
		q1_ = out2.q1;
		z0_ = out2.z0;
		k_ = [];
		[x__ y__ cflag__, k__] = bvp2wavetwopass(@out2.odefun,@out2.bcfun,Xi,opt);
	case {'wtz'}
		[x_ z1_ q1_ ] = wavetrainz(@cfun,Xi,omega,opt,hfun);
		z1_=z1_(:,1);
		q1_=q1_(:,1);
	case {'wtq'}
		out2 = out();
		% reset
		out2.Q1 = [];
		out2.solver = @bvp2wavetrain;
		out2.solve();
		x_  = out2.x;
		z1_ = out2.z1;
		q1_ = out2.q1;
		z0_ = out2.z0;
%		[x_ z1_ q1_ ] = wavetrainq(@cfun,Xi,az1,omega,hfun,wfun);
%		z1_=z1_(:,1);
%		q1_=q1_(:,1);
	otherwise
		error('here');
	end % switch
	if (isempty(k_))
		k_ = 1./q1_.*cdiff(q1_,[],2)./cdiff(x_,[],2);
	end
	if (0)
		scale = 1/z1_(1);
	else
		scale = 1/az1;
	end
	z1(:,jdx) = scale*interp1(x_,z1_,x);
	q1(:,jdx) = scale*interp1(x_,q1_,x);
	z0(:,jdx) = scale*interp1(x_,z0_,x);
	k(:,jdx) = interp1(x_,k_,x);
	%k__ = interp1(x_,k_,x);
    end % for jdx

%r = damping_modulus_river(Q0,w,h0,cd,omega);
%[r_ k_] = damping_modulus(Q0,w0,h0,cdfun(0),omega,az1);
%[r_(2) k_(2)] = damping_modulus(Q0,wfun(Xi(2)),hfun(Xi(2)),cdfun(Xi(2)),omega,0); %az1);
c = out.odefun(out.x,out.Q1);
R = roots2(c);
r_ = imag(R(:,1));
k_ = real(R(:,1));
%figure()
%clf
%plot([imag(R),real(R)])
%pause
	

figure(idx);
clf
subplot(3,2,1)
plot(x,abs(z1),'linewidth',1);
hold on;
if (1==idx)
	plot(x(1:1:end),abs(exp(-r_*x(1:1:end))),'.');
end
ylim([0 2])
ylabel('|z_1|');

subplot(3,2,2);
plot(x,wrapTo2Pi(angle(z1)));
hold on
ylabel('arg(z_1)');
%plot(x,wrapTo2Pi(k_*x),'.');
legend(leg{:})

subplot(3,2,3);
plot(x,abs(q1),'linewidth',1);
hold on;
%plot(x(1:1:end),abs(exp(-r_*x(1:1:end))),'.');
ylim([0 sqrt(g*h0)])
ylabel('|q_1|');

subplot(3,2,4);
plot(x,wrapTo2Pi(angle(q1)));
ylabel('arg(q_1)');
hold on
if (1==idx)
plot(x,wrapTo2Pi(k_*x),'.');
end

subplot(3,2,5);
u = bsxfun(@times,q1,1./hfun(x));
plot(x,abs(u),'linewidth',1);

figure(10+idx);
clf();
subplot(2,2,1)
if (0)
c=cfun(x(1:end/2));
Km = max(abs(sqrt(c(:,3))));
end
plot(x,[-real(k) -real(R(:,1) -real(k__))]);
%hline(r_,'g')
%ylim(Km*[-1 1])

subplot(2,2,2)
plot(x,[imag(k) imag(R(:,1)) imag(k__)]);
%hline(k_,'g')
k0 = omega/sqrt(g*h0);
hline(k0,'b');
%ylim(Km*[-1 1])

subplot(2,2,3)
h=hfun(x);
q1(1)
h(1)
%plot(x,[(h./g).*(abs(q1(:,1))./h).^2, g*abs(z1(:,1)).^2]);
%plot(x,[(abs(q1(:,1))./h)*NaN (abs(z1(:,1)).*sq1rt(g./h)).^2]);
lambda = 2*pi*sqrt(g.*hfun(x))./omega;
E = [h.*(abs(q1(:,1))./h).^2 g*(abs(z1(:,1))).^2];
plot(x,[lambda.*sum(E,2),lambda.*E(:,1), lambda.*E(:,2)])
%plot(x,[(1./h.*abs(q1(:,1))).^2/g g*(abs(z1(:,1))).^2]);
%plot(x,[(abs(q1(:,1))./h).^2/g g*(abs(z1(:,1))).^2]);
%plot(x,[(h.*abs(q1(:,1))./h).^2 + abs(z1(:,1)).^2]);
%plot(x,[0.*(abs(q1(:,1))./h).^2 + abs(z1(:,1)).^2]);
title('energy')

figure(100);
subplot(2,2,idx);
plot(x,hfun(x));
ylabel('h0')
subplot(2,2,idx+2);
plot(x,cdfun(x));
ylabel('cd')

figure(101)
subplot(2,2,1)
plot(x,z0)
ylabel('z0');
subplot(2,2,2)
plot(x,wfun(x))
ylabel('w');
subplot(2,2,3)
plot(x,abs(out.Q1))
ylabel('Q1');
subplot(2,2,4)
plot(x,abs(out.q1))
ylabel('q1')

figure(1001);
%plot(x,abs([ ((z1(:,1)-0*z1(1))/abs(z1(1))), -(q1(:,1)-0*q1(1))/abs(q1(1))]))
%plot(x,[1/(p-1)*(abs(z1(:,1))/abs(z1(1))-1), 1./(1-1/p)*(1-abs(q1(:,1))/abs(q1(1)))])
plot(x,[ abs(z1(:,1))/abs(z1(1)), abs(q1(1))./abs(q1(:,1))])
%ylim([1 p.^(1.1)]);
ylim([0 2.1]);

figure(1002);
c  = out.odefun(out.x,out.Q1);
c_ = c(:,3)./c(:,1);
c  = cfun(x,out.Q1);
subplot(2,2,1)
plot(real([c_ c(:,3)]))
subplot(2,2,2)
plot(imag([c_ c(:,3)]))

end % for idx

% case 2: depth convergence w/o daping

function c  = cfun(x,Q1)
	h_  = hfun(x);
	cd_ = cdfun(x);
	w   = wfun(x);
	dx  = L/(opt.nx-1);
	dw_dx = 1/dx*(wfun(x+0.5*dx)-wfun(x-0.5*dx));
	q0  = Q0./w;
	qhr = abs(Q1./w);
	cf  = friction_coefficient_dronkers(cvec(q0./qhr),2);
	c1  = cf(:,2);
	c2  = cf(:,3);
	c   = ones(length(x),3);
	c(:,2) = -1./w.*dw_dx;
	c(:,3) = ( omega.^2./(g.*h_) ...
	      - 1i*omega*(2*cd_.*c2.*q0)./(g*h_.^3*pi) ...	% river damping
	      - 1i*omega*(cd_.*c1.*qt0)./(g*h_.^3*pi) ...	% self damping
	    );
end

function y = fun(x,val0,valL)
%	y = val0 + (valL-val0).*(x-0.5*L)/(0.5*L).*(x>0.5*L);
	% quadratic
	y = val0 + (valL-val0).*((x-0.5*L)/(0.5*L)).^2.*(x>0.5*L);
end

end % test_wave_twopass


