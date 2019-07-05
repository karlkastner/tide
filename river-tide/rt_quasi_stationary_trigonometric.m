% Sat 30 Sep 14:54:34 CEST 2017
%% quasi statinary form of the SWE
function [out] = swe_quasi_stationary(X,w,Q0,cd,omega,fzb,bc,opt)
	nx = opt.nx;
	nf = opt.nf;

	q0 = Q0/w;
	sq0 = sign(q0);

	g = Constant.gravity();

	% discretise space
	dx = (X(2)-X(1))/(nx-1);
	x  = X(1) + dx*(0:nx-1)';

	% bed level
	zb = fzb(x);

	if (isfield(opt,'aa'))
		aa = opt.aa;
	else
		aa = 1;
	end

	% solver options
	if (isfield(opt,'maxiter'))
		sopt.maxiter = opt.maxiter;
	else
		sopt.maxiter = 100;
	end
	% relaxation
	if (isfield(opt,'relaxation'))
		sopt.p = opt.relaxation;
	else
		sopt.p       = 0.01;
	end
	
	sopt.verbose = 1;
	sopt.reltol  = 1e-12;

	% initial values
	%hupstream   = normal_flow_depth(Q0,W0,50,S);
	% TODO choose normal flow depth and linearly interpolate z1
	for idx=0:length(nf)
	switch (idx)
	case {0}
		z0 = max(0,zb+20);
		y  = z0;
	case {1}
		cz1 = 0*rand(nx,1);
		sz1 = 0*rand(nx,1);
		cq1 = 0*rand(nx,1);
		sq1 = 0*rand(nx,1);
		% stack
		y = [y; cz1; sz1; cq1; sq1];
	otherwise {2}
		error ('not yet implemented');
	end % switch
	end %for

	% solve non-linear system by picard iteration
	y = picard(@fun,y,sopt);

	% extract and write
	out = extract(y,nx,nf,true,true);
	out.x  = x;
	out.zb = zb;
	out.h0 = out.z0-zb;
	out.opt = opt;

	out.w     = w;
	out.Q0    = Q0;
	out.cd    = cd;
	out.omega = omega;

	function y = fun(y)
		
		% extract
		out = extract(y,nx, nf,true);
%		figure(1)
%		subplot(2,2,1)
%		plot(x,[zb out.z0])
%		subplot(2,2,2)
%		plot(x,out.az1)
%		subplot(2,2,3)
%		plot(x,out.aq1)
%pause(0.1)
		if (1)
			in = [1:nx-1];
		else
			in = [2:nx];
		end
		in = 2:nx-1;
		sc = @(x) select(x,in,2);

		z0  = diag(sparse(out.z0));
		cz1 = diag(sparse(out.cz1));
		sz1 = diag(sparse(out.sz1));
		cq1 = diag(sparse(out.cq1));
		sq1 = diag(sparse(out.sq1));

		h0 = out.z0-zb;
%min(h0(2:end))
		ih0 = diag(sparse(1./h0));
		h0  = diag(sparse(h0));

		z0   =  z0(in,:);
		cz1  = cz1(in,:);
		sz1  = sz1(in,:);
		cq1  = cq1(in,:);
		sq1  = sq1(in,:);
		h0   =  h0(in,:);
		ih0  =  ih0(in,:);

		if (0)
			alpha = abs(q0)./abs(q1);                                               
        	        p     = -friction_coefficient_dronkers(alpha);
		end

		% central differences for the n-1 interior points
		Dc = 0.5/dx*spdiags(ones(nx,1)*[-1 0 1],-1:1,nx,nx);
		Dc = Dc(in,:);
		%Dci = in(Dc);

		% forward differences
		Dr  = 1/dx*spdiags(ones(nx,1)*[-1 1],0:1,nx,nx);
		Dr = Dr(in,:);
%		Di = in(D);
%		Dr = Dr(2:end,:);

		Dl = 1/dx*spdiags(ones(nx,1)*[-1 1],-1:0,nx,nx);
		Dl = Dl(in,:);
%		Dl = Dl(2:end,:);
%		Dl = Dr;
%		D = Dc;
%		Dl = Dc;

		if (sq0>0)
			D = Dr;
		else
			D = Dl;
		end
%		D=Dc;

		I  = eye(nx,nx);
		I  = I(in,:);
		%II = speye(nx);
		%Z  = zeros(nx-2,nx);

% construct discretisation matrix

		z = zeros(length(in),1);
		Z = zeros(length(in),nx);
	
		% continuity c1
		Acc1 = [ Z, ... % z0
			 Z, ... % cz1
			 omega*I, ... % sz1
		         Dc, ... % cq1
			 Z      % sq1
			];
		bcc1 = z;
	
		% continuity s1
		Acs1 = [ Z, ... % z0
			-omega*I, ... % zc1
			 Z, ... % zs1
			 Z, ... % cq1
			 Dc      % sq1
			];
		bcs1 = z;
	
		% zero order term
		Am0 = [	  sc(g*h0 - (aa*q0.^2).*ih0.^2 - (aa*sq1.^2.*ih0.^2/2) - (aa*cq1.^2.*ih0.^2/2))*D, ...
	  		+ sc(1/2*cz1*g - aa*q0*cq1.*ih0.^2)*Dc, ...
			+ sc(-aa*q0*sq1.*ih0.^2 + 1/2*g*sz1)*Dc, ...
			+ (cd*sq0*cq1).*ih0.^2/2 + sc(aa*cq1.*ih0)*Dc, ...
		        + (cd*sq0*sq1).*ih0.^2/2 + sc(aa*sq1.*ih0)*Dc
		      ];
		bm0 = -(diag(sc((cd*q0.^2*sq0).*ih0.^2)) + sc(aa*sq1.^2.*ih0.^2/2 + (aa*cq1.^2.*ih0.^2)/2 + (aa*q0.^2.*ih0.^2))*D*zb);
	
		% first order cosine
		Amc1 = [  sc(cz1*g - 2*aa*cq1*q0.*ih0.^2)*Dc, ...
			+ sc(g*h0 - (aa*q0.^2.*ih0.^2) - (aa*sq1.^2.*ih0.^2)/4 - (3*aa*cq1.^2.*ih0.^2)/4)*Dc, ...
			- sc(aa*cq1.*sq1.*ih0.^2)*Dc/2 , ...
 			+ (2*cd*q0*sq0).*ih0.^2 + sc(2*aa*q0.*ih0)*Dc, ...
			+ omega*I
			];
		bmc1 = -sc(2*aa*cq1*q0.*ih0.^2)*Dc*zb;
	
		% first order sine
		Ams1 = [  sc(g*sz1 - 2*aa*q0*sq1.*ih0.^2)*Dc, ...
			- sc(aa*cq1.*sq1.*ih0.^2)*Dc/2 , ...
			+ sc(g*h0 - (aa*q0.^2.*ih0.^2) - (3*aa*sq1.^2.*ih0.^2)/4 - (aa*cq1.^2.*ih0.^2)/4)*Dc, ...
			- omega*I, ...
			+ (2*cd*q0*sq0).*ih0.^2 + sc(2*aa*q0.*ih0)*Dc 
		];
		bms1 = -sc(2*aa*q0*sq1.*ih0.^2)*Dc*zb;
	
		%
		% downstream boundary condition
		%
		z = zeros(1,nx);
		A_bc = [[1,z(1:nx-1)],            z,            z,                    z,                   z; %  z0
	                           z, [1,z(1:nx-1)],            z,                    z,                   z; % cz1
			           z,            z, [1,z(1:nx-1)],                    z,                   z; % sz1
			           z,            z,            z, [[-1,1]/dx,z(1:nx-2)],                   z; % cq1
			           z,            z,            z,                    z, [[-1,1]/dx,z(1:nx-2)] % sq1
			];
		b_bc = [ bc.z0(1);
		         bc.cz1(1);
		         bc.sz1(1);
		        -omega*bc.sz1(1);
		        +omega*bc.cz1(1)];
		%
		% upstream boundary condition
		%
		% (no mean flow bc necessary upstream)
		A_bc = [A_bc;
			[z(1:nx-1),1],            z,            z,                    z,                   z; % z0
	                           z, [z(1:nx-1),1],            z,                    z,                   z; % cz1
			           z,            z, [z(1:nx-1),1],                    z,                   z; % sz1
			           z,            z,            z, [z(1:nx-2),[-1,1]/dx],                   z; % cq1
			           z,            z,            z,                    z, [z(1:nx-2),[-1,1]/dx] % sq1
			];
%bc.z0
%out.z0(end-1:end)
%ause
		b_bc = [ b_bc;
			 bc.z0(2);
		         bc.cz1(2);
		         bc.sz1(2);
		        -omega*bc.sz1(2);
		        +omega*bc.cz1(2)];
	%	A_bc = [A_bc;
	%		z,z 1, z z z;
	%		z,z,z 1, z z;
	
		% stack
		A = [Acc1;
	             Acs1;
	              Am0;
	             Amc1;
	             Ams1
		     A_bc];
	
		b = [bcc1;
	             bcs1;
	             bm0;
	             bmc1;
	             bms1;
		     b_bc];

%figure(11)
%clf
%	spy(A)
%	spy(Am0)
%pause

		% solve
		y = A \ b;

		% limit bed elevation
		z0 = y(1:nx);
		hmin=1e-2;
		z0 = max(z0,zb+hmin);
%		z0 = max(zb+1e-3*(out.z0-zb),z0);
		y(1:nx) = z0;
		
		% limit
%		y(1:nx) = max(zb+hmin,y(1:nx));
	end % fun

	function out = extract(y,nx,nf,aflag,pflag)
		out = struct();
		for idx=0:nf
		switch (idx)
		case {0}
			out.z0  = y(1:nx);
		case {1}
			out.cz1 = y(nx+1:2*nx);
			out.sz1 = y(2*nx+1:3*nx);
			out.cq1 = y(3*nx+1:4*nx);
			out.sq1 = y(4*nx+1:5*nx);
		if (nargin() > 3 && aflag)
			out.az1 = hypot(out.cz1,out.sz1);
			out.aq1 = hypot(out.cq1,out.sq1);
		end
		if (nargin() > 4 && pflag)
			out.pz1 = atan2(out.sz1,out.cz1);
			out.pq1 = atan2(out.sq1,out.cq1);
		end
		otherwise
			error('not yet implmented');
		end % switch
		end % for
	end % function extract
end % swe_quasi stationary


