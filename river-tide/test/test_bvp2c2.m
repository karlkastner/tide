% Wed 13 Mar 15:36:35 CET 2019

% bvp2c2 is second order accurate, bvp2c only first order accurate
%
% bvp2c2 is only secon order accurate, if the water depth changes smoothly along
% the channel, or if the grid points coincide with steps (!)
%
%
function [x,y1,y2,ypm5,A1,A2,b1,b2] = test_bvp2c2
	w     = 500;
	omega = 2*pi/86400;
	cd    = 0.0025;
% NOTE: Q0<0, as discharge is directed into the downstream direction
	Q0 = -1e3;
	L  = 3e5; %4e5;
	X  = [0, L];
	g  = 9.81;
	h0 = 15;
	%dh_dx = -(h0-hu)/(L2);
	L2     = L; %2.9e5;
	%L2     = 0.5*L;
	dzb_dx = 15/L2;
	hu     = normal_flow_depth(Q0,w,sqrt(g/cd),dzb_dx);
%	cd  = 0.01*cd;
%	hu     = 5; %0.1;
	hu     = 0.1*h0;
	dh_dx  = (-dzb_dx + hu/L2);
	%hu     = %0.5*h0;
	hfun   = @(x) min(h0,max(hu,h0 + dh_dx*x));
	%dh_dx = 0;
	%hfun   = @(x) hu+(h0-hu)*((L-x)/L).^2;
	%hfun = @(x) h0*exp(-x/L);

% -> TODO compute discharge and from there water level

% plotting
if (0)

	disp('characteristic polynomial coefficients')
	c = odefun(0);
	%disp('c(1:3))
	%num2str(c)
	fprintf('%+3.1e %+3.1ei\n',[real(c(:)),imag(c(:))]')
	disp('characteristic polynomial roots')
	roots(c(1:3))
%pause

	bc        = struct();
	bc(1).val =  1;
	bc(1).p   = [1,0;
                     0,0];
	bc(2).val = 0;
	bc(2).p   = [0,0;
                     1,0];

	opt    = struct();
	opt.nx = 400;
	%opt.sopt.miniter = 100; %opt.nx;
	opt.sopt.reltol = 1e-12;
	xr     = linspace(0,L)';
	opt.xr = xr;

	figure(1);
	clf();
	figure(2);
	clf();
	figure(3);
	clf();

	id = 0;
	leg={};

	%[x, y, cflag, dydx, l, c, A]
	id = id+1;
	[ x{id},  y{id}, ypm{id}, v1, v2, v3, v4, A1, b1, yr] = bvp2c( @odefun, @bcfun, X, opt);
	ypm{id} = [ypm{id}(:,1),ypm{id}(:,3)];
	xc{id}  = mid(x{id});
	leg{id} = 'bvp2c'; 
if (0)	
	id = id+1;
	[x{id}, y{id}, ypm{id}, A2, b2,yr] = bvp2c2(@odefun,    bc, X, opt);
	ypm{id} = [ypm{id}(:,2),ypm{id}(:,1)];
	xc{id} = mid(x{id});
	leg{id} = 'bvp2c2'; 
end
if (0)
	% solve ypm as ode
	id = id+1;
	opt_ = odeset('maxstep',0.01*L/opt.nx,'reltol',sqrt(eps),'abstol',sqrt(eps));
	%y0 = 1.5*ypm1(1,:) - 0.5*ypm1(2,:);
	y0    = ypm{1}(1,:);
	[xc{id},ypm{id}] = ode23s(@dQlr_dx,[mid(x{1}(1:2)),L],y0([1,2]),opt_);
	y{id} = sum(ypm{id},2);
	x{id} = xc{id};
	leg{id} = 'ode lr'; 

	% solve y ode, with Q and dQ/dx
	id = id+1;
	% accuracy is paramount here to avoid wiggles
	d1 = (y{1}(2)-y{1}(1))./(x{1}(2)-x{1}(1));
	d2 = (y{1}(3)-y{1}(2))./(x{1}(3)-x{1}(2));
	y0 = [y{1}(1); 1.5*d1 - 0.5*d2];
	%(y{1}(2)-y{1}(1))./(x{1}(2)-x{1}(1))];
	[x{id},y{id}] = ode23s(@dQ_dx,[0,L],y0,opt_);
	y{id} = y{id}(:,1);
	ypm{id} = NaN(length(y{id}),2);
	xc{id} = x{id};
	leg{id} = 'ode dot'; 
end

if (0)
	constflag = [false,true];
	% for idx=1:2
	for idx=1:2
	id = id+1;
	[y{id}, ypm{id}] = rt_simple(Q0,h0,dzb_dx,x{1},constflag(idx));
	x{id}  = x{1};
	xc{id} = x{id};
	leg{id} = 'simple'; 
	end
end
	id = id+1;
	[x{id},y{id}, xc{id}, ypm{id},dy_dx{id}] = rt_simple2(@odefun,x{1});
	leg{id} = 'simple';

%$b1(2:3:end) = [];
%[y1, y2]
%[y1./y2]
%[b1,b2]
%A1(2:3:end,:)=[];A1(:,2:3:end)=[]; A1=swapevenodd(A1.').'; norm(A1-A2) %full(A1),full(A), norm(full(A1-A))
	
	for idx=1:length(x)
		figure(1);
		subplot(2,2,1)
		plot(x{idx},abs(y{idx}));
		title('Amplitude');
		hold on
		if (length(x)==idx)
			legend(leg{:}) %'bvp2c','bvp2c2','odeLR','ode','approx');
		end
		ylim([0,2]);

		subplot(2,2,2)
		plot(x{idx},angle(y{idx}));
		hold on
		title('Phase')

		figure(2);
		subplot(2,2,1);
		plot(xc{idx},abs(ypm{idx}(:,1)));
		hold on;
		title('Abs Right going');
		subplot(2,2,2);
		plot(xc{idx},abs(ypm{idx}(:,2)),'linewidth',2);
		hold on;
		title('Abs Left going');
		subplot(2,2,3);
		plot(xc{idx},angle(ypm{idx}(:,1)));
		hold on;
		title('Angle Right going');
		subplot(2,2,4);
		plot(xc{idx},angle(ypm{idx}(:,2)),'linewidth',2);
		hold on;
		title('Angle Left going');

		figure(3);
		subplot(2,2,1)
		plot(xc{idx},real(ypm{idx}(:,1)));
		hold on
		title('Right going real');
		subplot(2,2,2)
		plot(xc{idx},real(ypm{idx}(:,2)));
		hold on
		title('Left going real');
		subplot(2,2,3)
		plot(xc{idx},imag(ypm{idx}(:,1)));
		hold on
		subplot(2,2,4)
		plot(xc{idx},imag(ypm{idx}(:,2)));
		hold on
	end

	figure(4)
	clf
	plot(x{1},-hfun(x{1}));

	figure(5);
	clf
	plot(x{end},[cdiff(y{end})./cdiff(x{end}),dy_dx{end}]);

	namedfigure(13,'Relative change of wave number');
	clf();
	c = odefun(x{1});
	c = roots2(c);
	subplot(2,2,1)
	plot(x{1},abs(c(:,1:2)))
	subplot(2,2,2);
	plot(x{1},abs(bsxfun(@times,bsxfun(@times,cdiff(c),1./cdiff(x{1})),1./(c(:,2)-c(:,1)))));
	subplot(2,2,3)
	plot(x{1},real(c(:,1:2)))
	subplot(2,2,4)
	plot(x{1},imag(c(:,1:2)))

else % if not plotting but convergence test

n0 = 0;
n=9;
opt.nx = 2^n+1;
N = (2.^(n0:n-1)).'+1;
%N = cvec(2:10);
[x_1, y_1, ypm, v1,v2,v3,v4,A1_, b1] = bvp2c(@odefun, @bcfun, X, opt);
%[x_2, y_2, ypm_, A2_, b2] = bvp2c2(@odefun, bc, X, opt);
[x_2,y_2] = rt_simple2(@odefun,x_1);
err = [];

namedfigure(3,'bvp2c')
clf
plot(x_1,abs(y_1),'.-');
hold on

namedfigure(4,'bvp2c 2')
clf
plot(x_2,abs(y_2),'.-');
hold on
x_1_   = x_1;
y_1_   = y_1;
x_2_   = x_2;
y_2_   = y_2;
opt.xr = x_1;

namedfigure(5,'bvp2c expanded');
clf
plot(x_1,abs(y_1));
hold on
%[y_1, y_2]
%b1
%b2
%pause
for idx=length(N):-1:1
 	%N(idx) = 2.^idx+1;
	opt.nx = N(idx);

	opt.sopt.miniter = opt.nx;

	[x1, y1, ypm, v1, v2, v3, v4, v5, v6, yr_1 ] = bvp2c(@odefun, @bcfun, X, opt);
	%err(idx,4) = norm(y_1_-yr_1)/sqrt(length(yr_1));
	if (1)
		y_1 = y_1(1:2:end);
		y_2 = y_2(1:2:end);
	else
		y_1 = interp1(x_1_,y_1_,x1,'spline');
		y_2 = interp1(x_2_,y_2_,x1,'spline');
	end

	figure(5);
	% plot(x1,abs(yr_1));

	figure(3);
	plot(x1,abs(y1),'.-');

	if (0)
		[x2, y2]       = bvp2c2(@odefun, bc, X, opt);
	else
		[x2,y2,xc,ypm] = rt_simple2(@odefun,x1);
	end

	err(idx,1) = norm(y1 - y_1)/sqrt(length(y1));

	err(idx,2) = norm(y2 - y_2)/sqrt(length(y2));
		
	% difference between solutions
	err(idx,3) = norm(y1-y2)./sqrt(length(y1));
	%sqrt(norm(y_1).^2+norm(y_2).^2);

	figure(4);
	plot(x2,abs(y2),'.-');
end

	figure(2);
	clf
	loglog(N,err,'.-');
	disp('convergence rate');
	A = [ones(size(N)),log(N)];
	c = A \ log(err);
	c
	hold on
	%loglog(N,exp(A*c))
	legend('bvp2c 1','bvp2c 2','difference'); 
end

% simplified river tide, assumes Qr > Qt
function c = odefun(x,~)
	if (0 == nargin())
		c = zeros(0,3,1);
	else
		nx = length(x);
		h  = hfun(x);
		%Area  = w*h;
		%c  = [1i*g*h/(omega*w),zeros(nx,1),cd./Area.^2*(1i*omega + 2*Q0),zeros(nx,1)];
		c  = [1i*g/(omega)*ones(nx,1),zeros(nx,1),(1i*omega./h + 2*cd*Q0./(w*h.^3)),zeros(nx,1)];
	end
end

function dQ_dx = dQ_dx(x,Q)
	c = odefun(x,Q);
	% c1 Q'' + c2 Q' + c3 Q = 0
	% Q' = Q2
	% c1 Q2' + c2 Q2 + c3 Q = 0
	A = [0, 1;
             -c(3)/c(1),-c(2)/c(1)];
	dQ_dx = A*Q;
end

% solve river tide as initial value problem
function dQlr_dx = dQlr_dx(x,Qlr)
	A = ode2_matrix(@odefun,x,Qlr,sqrt(eps)*L,'t');
	dQlr_dx = A*Qlr(:);
end


% solve river tide as nearly constant coefficient boundary value problem
% separate domain in two reaches, backwater affected and uniform flow reach
% propagate along each reach individually
function [y, ypm] = rt_simple(~,h0,dzb_dx,x,constflag)
	% TODO, this is not the true linearized backwater eq
	L       = h0/dzb_dx;
	xc_     = L/2;
	xc_     = 0;
	xc_     = L/2;
	dx_     = sqrt(eps)*L; % L/2
	[A, Ac] = ode2_matrix(@odefun,xc_,[],dx_,'t');
	if (nargin() > 4 && constflag)
		A = Ac;
	end
	[V,E] = eig(A);
	E     = diag(E);
%	E(1)=conj(E(1));
	%E(1)  = -E(1);
	%E     = [E(1);-E(2)];
	Vi    = inv(V);
	A0    = V*diag(exp((0-xc_)*E))*Vi;
	AL    = V*diag(exp((L-xc_)*E))*Vi;
	% unity incoming from right and zero incoming from left
	rhs   = [1;0];
	A     = [A0(2,:);
                 AL(1,:)];
	%A     = [y0(:,1),yr(:,2)];
	yhat  = A \ rhs;
%	A0*yhat
%	AL*yhat
%pause
	for jdx=1:length(x)
		% TODO simplify
		ypm(jdx,:) = V*diag(exp((x(jdx)-xc_)*E))*Vi*yhat;
%pause
	end
	if (~constflag)
	yL         = AL*yhat;
	[~, Au] = ode2_matrix(@odefun,L,[],[],'t');
%Au
%yL
%pause
	[V,E] = eig(Au);
	E     = diag(E);
	fdx=find(x>L);
	if (~isempty(fdx))
	for jdx=rvec(fdx)
		ypm(jdx,:) = (V*diag(exp(E*rvec(x(jdx)-L)))*yL);
	end
	end
%	ypm(x>L,:) = NaN;
	end

	y = sum(ypm,2);
end % rt_simple


% TODO switch by index
function [v,p,q] = bcfun(x,~,~)
	if (x == X(1))
		v = 1;
		p = [1,0];
		q = [0,1];
	else
		v = 0;
		p = [1,0];
		q = [1,0];
	end
end

end % test_bvp2c

%
% exterior functions
%

function [x,y,xc,ypm,dy_dx] = rt_simple2(odefun,x)
	mode.step = 'analytic';
	flag = true;
	%flag = false;
	%mode.step = 'fdm';
	ns  = length(x)-1;
	dx  = diff(x);
	xc  = mid(x);
	I   = eye(2);
	A   = zeros(2*ns);
	rhs = zeros(2*ns,1);

	% get system matrix
	D  = [];
	Dl = [];
	Dr = [];
	Al = [];
	Ar = [];
	for idx=1:ns
		[D,Dc] = ode2_matrix(odefun, xc(idx),[],dx(idx),'t');
		if (~flag)
			D = Dc;
		end
		switch (mode.step)
		case {'fdm'}
			% V exp(1/2*dx*E)V' ~ (1 + dx/2*D)
			% left and right wave at left end
			Al(:,:,idx) = (I - 0.5*dx(idx)*D); %(:,:,idx));
			% left and right wave at right end
			Ar(:,:,idx) = (I + 0.5*dx(idx)*D); %(:,:,idx));
			Dl(:,:,idx) = D; %(:,:,idx);
			Dr(:,:,idx) = D;% (:,:,idx);
		case {'analytic'}
			[V,E] = eig(D); %(:,:,idx));
			E     = diag(E);
			%E = -[E(1);-E(2)]; 
			%E = -conj(E);
			Vi    = inv(V);
			Al(:,:,idx)    = V*diag(exp(-0.5*dx(idx)*E))*Vi;
			Ar(:,:,idx)    = V*diag(exp(+0.5*dx(idx)*E))*Vi;
			% TODO why does this only work when convergence is
			%      neglected for the derivative at the interface?
			%Dl(:,:,idx)    = V*diag(E.*exp(-0.5*dx(idx)*E))*Vi;
			%Dr(:,:,idx)    = V*diag(E.*exp(+0.5*dx(idx)*E))*Vi;
			Dl(:,:,idx)    = V*diag(exp(-0.5*dx(idx)*E))*Vi*Dc;
			Dr(:,:,idx)    = V*diag(exp(+0.5*dx(idx)*E))*Vi*Dc;
		end
	end % switch

	neq = 1;
	% incoming wave
	A(neq,1) = Al(2,1,1);
	A(neq,2) = Al(2,2,1);
	% TODO user defined
	rhs(neq)   = 1; 

	% % for interior points, continuity
	% for each segment, value at left, value at right
	for idx=1:ns
		if (true) %false) %flag)
			% the left and right going wave cannot be contiuous between segments
			% when the change of the wave number is not considered,
			% as part of the wave is reflected at the interface,
			% but the total wave is continuous and has a continuous derivative
			%
			if (idx > 1)
				% continuity at left end, right going wave
				A(neq-1,2*idx-1) = -Al(1,1,idx);
				A(neq-1,2*idx)   = -Al(1,2,idx);
				% continuity at left end, left going wave
				A(neq,2*idx-1)   = -Al(2,1,idx);
				A(neq,2*idx)     = -Al(2,2,idx);
			end
			if (idx < ns)
				% continuity at right end, right going wave
				A(neq+1,2*idx-1)    = Ar(1,1,idx);
				A(neq+1,2*idx)      = Ar(1,2,idx);
				% continuity at right end, left going wave
				A(neq+2,2*idx-1)    = Ar(2,1,idx);
				A(neq+2,2*idx)      = Ar(2,2,idx);
				neq = neq+2;
			end		
		else
			if (idx>1)
				% continuity of value
				A(neq-1,2*idx-1)    = -(Al(1,1,idx) + Al(2,1,idx));
				A(neq-1,2*idx)      = -(Al(1,2,idx) + Al(2,2,idx));
				% continuity of derivative
				A(  neq,2*idx-1)    = -(Dl(1,1,idx) + Dl(2,1,idx));
				A(  neq,2*idx)      = -(Dl(1,2,idx) + Dl(2,2,idx));
			end
			if (idx < ns)
				% continuity of value
				A(neq+1,2*idx-1)    = (Ar(1,1,idx) + Ar(2,1,idx));
				A(neq+1,2*idx)      = (Ar(1,2,idx) + Ar(2,2,idx));
				% continuity of derivative
				A(neq+2,2*idx-1)    = (Dr(1,1,idx) + Dr(2,1,idx));
				A(neq+2,2*idx)      = (Dr(1,2,idx) + Dr(2,2,idx));
				neq = neq+2;
			end
		end
	end % for idx
	% for last point, right condition
	neq = neq+1;
	A(neq,2*ns-1) = Ar(1,1,ns);
	A(neq,2*ns)   = Ar(1,2,ns);
	% no incoming wave from right
	% TODO, user defined
	rhs(2*ns) = 0;

%	s  = 1./abs(diag(A));
%	[min(s),max(s)]
if (0)
	s = 1./max(abs(A)').';
	A = diag(s)*A;
	rhs = s.*rhs;	
end
%s
%full(abs(A))
%Ar(:,:,end)
%	[neq,2*ns-1]
%pause

	% solve
	ypm  = A \ rhs;
	% reshape
	ypm = reshape(ypm,2,[]).';
	%x = xc;
	y   = zeros(size(x));
	dy_dx = zeros(size(x));
	for idx=1:ns
		y(idx)     = sum(Al(:,:,idx)*ypm(idx,:).');
		dy_dx(idx) = sum(Dl(:,:,idx)*ypm(idx,:).');
		if (idx<ns)
		dy_dx(idx+1,2) = sum(Dr(:,:,idx)*ypm(idx,:).');
		end
	end
	y(ns+1) = sum(Ar(:,:,ns)*ypm(ns,:).');
	dy_dx(ns+1) = sum(Dr(:,:,ns)*ypm(ns,:).');
	%y   = sum(ypm,2); 
	%y = inner2outer(y);
end % rt_simple2

