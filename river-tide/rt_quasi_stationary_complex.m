% Fri  8 Sep 08:59:04 CEST 2017
%% quasi-stationary solution of the SWE
%% TODO staggered grid does not help: q1' needed
% 
function [out] = rt_quasi_stationary(X,w,Q0,cd,omega,fzb,bc,opt)
	nx = opt.nx;
	nf = opt.nf;
	% advective acceleration flag
	f.aa = opt.flag.aa;

	q0 = Q0/w;
%	w=1;

	g = Constant.gravity();

	% discretise space
	dx = (X(2)-X(1))/(nx-1);
	x  = X(1) + dx*(0:nx-1)';
	% xc = x(1:end-1)+0.5*dx;

	% get bed level
	zb = fzb(x);

	sopt.maxiter = 200;
	% relaxation
	sopt.p = 0.125;
	sopt.verbose = 1;
	sopt.reltol = 1e-12;

	% initial values
	%hupstream   = normal_flow_depth(Q0,W0,50,S);
	% TODO choose normal flow depth and linearly interpolate z1
	z0 = max(0,zb+10); % + rand(nx,1);
	y = z0;

	if (opt.nf>0)
		z1 = 1e-7*(rand(nx,1)+1i*rand(nx,1));
		q1 = 1e-7*(rand(nx,1)+1i*rand(nx,1));
		% stack
		y = [y; z1; q1];
	end
%figure()
%clf
%plot([zb z0])
%pause

	% solve non-linear system by picard iteration
	y = picard(@fun,y,sopt);

	% extract and write
	out.x  = x;
	out.zb = zb;
	z0 = y(1:nx);
	out.z0 = z0;
	out.h0 = h0;
	out.q0 = q0;
	out.opt = opt;

	if (nf > 0)
		z1 = y(nx+1:2*nx);
		q1 = y(2*nx+1:3*nx);
		z1 = cmean(z1);
		q1 = cmean(q1);
		out.z1 = z1;
		out.q1 = q1;
	end

	function y = fun(y)
		
		% extract
		z0 = y(1:nx);
		if (nf > 0)
			z1 = y(nx+1:2*nx);
			q1 = 1e-7*y(2*nx+1:3*nx);
			aq1 = abs(q1);
		else
			q1 = 0;
			aq1 = sqrt(eps);
		end

%		aq1 = max(aq1,1e-7);

		% smooth
%		z0 = cmean(z0);

		h0  = z0-zb;

		alpha = abs(q0)./abs(q1);                                               
                p     = -friction_coefficient_dronkers(alpha);

%		figure()
%		clf()
%		plot([h0 abs(z1)])
%		hold on
%		pause(0.1)

		% central differences for the n-1 interior points
		Dc = 0.5/dx*spdiags(ones(nx,1)*[-1 0 1],-1:1,nx,nx);
		Dci = in(Dc);
		% forward differences, 
		D = 1/dx*spdiags(ones(nx,1)*[-1 1],0:1,nx,nx);
		Di = in(D);
		%Dc = 1/dx*spdiags(ones(nx,1)*[-1 1],-1:0,nx,nx);
		%Dci = in(Dc);
%		Dc = D; Dci=Di;
%		full(Di)
%pause
		%D = D(2:end-1,:);
		I = ones(nx,1);
		II = speye(nx);

		Z = zeros(nx-2,nx);

	
		% momentum 0:
		% (cd (aq1^2 p0 + aq1 p1 q0 + p2 q0^2 + (p3 q0^3)/aq1))/(h0^2 Pi)
		%  + g h0 D z0 
		%  - (q0^2 (D z0 - D zb))/h0^2

		% discretisation matrix
		A_m0 = [sd(in(g*h0 - f.aa*q0^2./h0.^2))*Di];
		% right hand side
		b_m0 = -in(cd./(pi*h0.^2).*(aq1.^2.*p(:,1) + aq1.*p(:,2)*q0 + p(:,3).*q0^2 ...
			  + p(:,4)*q0^3./aq1 )) ...
		          - f.aa*in(q0.^2./h0.^2).*(Di*zb);
		% bc for water surface elevation
		% dirichlet for ws, neuman for discharge
		A_upstream   = [ 1, zeros(1,1*nx-1) ];
		A_downstream = [zeros(1,nx-1),1 ];
		b_upstream   = [bc.upstream.z0 ];
		b_downstream = [bc.downstream.z0 ];

		if (nf>0)
			A_m0 = [A_m0,Z,Z];
			A_upstream = [A_upstream,zeros(1,2*nx)];
			A_downstream = [A_downstream,zeros(1,2*nx)];
		end

		A = A_m0;
		b = b_m0;

		
		% momentum 1
%		i o1 q1 
%		+ cd ((aq1 p1 q1)/(h0^2 pi) 
%		  + (2 p2 q0 q1)/(h0^2 pi) 
%		  + (3 p3 q0^2 q1)/(aq1 h0^2 pi)) 	
% i		- ( 2 cd ((aq1^2 p0)/(h0^2 pi) 
% i			+ (aq1 p1 q0)/(h0^2 pi) 
% i			+ (p2 q0^2)/( h0^2 pi) 
% i			+ (p3 q0^3)/(aq1 h0^2 pi)) z1/h0)
%		+ ( 2 D q1)/h0 
%		+ g z1 D z0 
%		+ g h0 D z1 
%		- (q0^2 D z1)/h0^2 
%		- ( 2 q0 q1 (D z0 - D zb))/h0^2 
%		+ ( 2 q0^2 z1 (D z0 - D zb))/h0^3 

% i o1 q1 
% + (cd (aq1 p1 + 2 p2 q0 + (3 p3 q0^2)/aq1) q1)/(h0^2 Pi) 
% - ( 2 cd (aq1^2 p0 + aq1 p1 q0 + p2 q0^2 + (p3 q0^3)/aq1) z1)/( h0^3 Pi)
% + (2 D q1)/h0 
% + g z1 D z0  		(!)
% + g h0 D z1 
% - (q0^2 D z1)/h0^2 
% - ( 2 q0 q1 (D z0 - D zb))/h0^2 	(!)
% + ( 2 q0^2 z1 (D z0 - D zb))/h0^3

		if (nf > 0)
			% continuity 1
			A_c1 = [Z, in(sd(1i*omega*I)), Dci];
			b_c1 = zeros(nx-2,1);

			% momentum 1
			A_m1 = [ [sd(in(g*z1))*Dci], ...
				 [sd(in(g*h0 - q0^2./h0.^2))*Dci ...
				 	+ in(sd((2*q0^2./h0.^3.*(Dc*(z0 - zb))))) ...
				  - in(sd(2*cd./(pi*h0.^3).*(p(:,1).*aq1.^2 + p(:,2)*q0 + p(:,3).*q0^2 + p(:,4).*q0^3./aq1))) ], ...
				 [ in(sd(1i*omega*I)) ...
					+ sd(in(2./h0))*Dci ...
					- in(sd((2*q0./h0.^2).*(Dc*(z0-zb)))) ...
				+ in(sd(cd./(pi*h0.^2).*(aq1.*p(:,2) + 2*p(:,3)*q0 + 3*p(:,4).*q0^2./aq1))) ]];
			b_m1 = zeros(nx-2,1);
			A_upstream = [A_upstream;
                        	         zeros(1,nx), 1, zeros(1,2*nx-1);
					 zeros(1,2*nx), -1/dx, 1/dx, zeros(1,nx-2)];
			A_downstream = [A_downstream;
	                                zeros(1,2*nx-1),1,zeros(1,nx);
					zeros(1,3*nx-2),-1/dx,1/dx];
			% rhs
			% dirichlet bc
			dz1_dt       = 1i*omega*bc.upstream.z1;
			b_upstream   = [b_upstream;
	       	                        bc.upstream.z1;
       		                        -dz1_dt];
			dz1_dt       = 1i*omega*bc.downstream.z1;
			b_downstream = [b_downstream;
	               	                bc.downstream.z1;
					-dz1_dt];
			A = [A; A_c1; A_m1];
			b = [b; b_c1; b_m1];
		end

		% stack
		A = [A; A_downstream; A_upstream];
		b = [b; b_downstream; b_upstream];
		
%		full(A)
%		full(b)
%pause		

		% solve
		y = A \ b;
		
		% limit
		hmin=1e-3;
		y(1:nx) = max(zb+hmin,y(1:nx));
		
	end % fun

	function X = in(X)
		X = X(2:end-1,:);
	end
	function sd = sd(X)
		sd = diag(sparse(X));
%		sd = sd(2:end-1,:);
	end
%	function sd = sd_(X)
%		sd = diag(sparse(cvec(X)));
%		sd = sd(2:end-1,2:end-1);
%	end
%	function sd = sd_(X)
%		sd = diag(sparse(cvec(X)));
%	end
end % swe_quasi stationary


