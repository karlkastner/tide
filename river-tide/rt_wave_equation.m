% Tue 11 Apr 13:38:18 CEST 2017
% Karl Kastner, Berlin
%
%% solve river tide as boundary value problem
%%
%% input:
%% omega : [nfx1] angluar frequency of tidal component, zero for mean flow
%% reach : [nrx1] struct
%% .L    : [1x1]  length of reaches
%%	 .width(x,h)	width
%%	 .bed(x,h)	bed level
%%	 .surface(x,h)	surface elevation
%%	 .Cd(x,h)	drag coefficient
%% .bc   : [nd,nf] boundary/junction conditions
%%	  bc(id,if).type : {surface, velocity, discharge} (dirichlet)
%%	  bc(id,if).val  : value
%% opt   : [1x1] struct
%%	- constant surface elevation
%%	- deactivative advective acceleration
%% 	.dx : spatial resolution
%%
%% dimensions:
%%	nr : nurmber or reaches
%%	nd : upstream/downstream index
%%	nf : frequency index
%
% TODO is a staggered grid necessary? -> good for bifurcation
% put elevation on grid point and add rows making [1,-1][z1,z2] = 0 and [1,-1],[z1,z3] = 0
% TODO omega must be integer multiple
function [x y] = river_tide(omega,reach,opt)
	% TODO no magic numbers
	z0 = 0;
	u0 = 1;
	if (nargin()<3)
		opt = struct();
	end
	if (~isfield(opt,'maxiter'))
		opt.maxiter = 10;
	end
	if (~isfield(opt,'abstol'))
		opt.abstol = 1e-3;
	end
	if (~isfield(opt,'dx'))
		opt.dx     = sum([reach.L])/100;
	end
	g = Constant.gravity();

	nf = length(omega);
	nr = length(reach);
	
	% preparation
	N = zeros(nr,1);
	for ir=1:nr
		% discretise space
		n = max(2,round(reach(ir).L/opt.dx));
		N(ir) = n;
		% discretise space
		Xr{ir} = reach(ir).L*(0:n-1)'/n;

		% TODO shortcut for constant functions
		% initial condition
		% TODO set non-zero depth
		z = z0*ones(n,1);
		u = u0*ones(n,1);
		y = [z; u];
	end % for ir

	iter = 0;
	% picard iteration
	while (1)
		% save last value (for convergence check)
		yold = y;

		% set up matrices for each reach
		for ir=1:nr
		n = N(ir);
		x = Xr{ir};
		
		% channel properties

		% bed level
		zb = reach(ir).bed(x);
			
		% width
		w = reach(ir).width(x);

		% drag coefficient
		cd = reach(ir).Cd(x);

		% for each frequency component
		for idf=1:nf
			% fetch unknowns
			% suface elevation
			% TODO reach local index
			zs = y(1:n);

			% velocity
			u  = y(n+1:2*n);

			% derived quantities

		 	% evaluate dynamic functions (friction, advective acceleration)
			% TODO dronkers
			r = 8/(3*pi)*u;

			% flow depth
			d = zs - zb;

			% area
			% assumes rectangular cross section
			A = w.*d;

			% differemce matrix
			D = derivative_matrix_1_1d(n,reach(ir).L);

			% TODO collect in buffers
			% TODO advective acceleration still missing
			% this is wrong, see derivation
			k=1;
			Mzz = [        ds(A)*D + ds(D*A)];
			Mzu = [     ds(1i*omega(idf)*k*w)];	% omega*k
			Muz = [ds(1i*omega(idf)*k + r./d)];	% omega*k
			Muu = [                  g*D];

			% right hand side
			rz = zeros(n,1);
			ru = zeros(n,1);

		 	% apply boundary condition
			% set dirichlet boundary conditions
			% TODO avoid deletion by setting up the matrix only for interior vertices
			apply_bc(1,reach(ir).bc(1,idf));
			apply_bc(n,reach(ir).bc(2,idf));

			% stack equations
			M = [Mzz, Mzu;
                             Muz, Muu];
			r = [rz; ru];
		end % for k (each frequency component)
		end % for jr (each reach)
		% TODO link reaches

		% solve system
		% TODO use sparse solver
		y = M \ r;
%		norm(y)
%		eigs(M,3,'LM')
		figure(iter+1)
		plot(x,[abs(y(1:n)) abs(y(n+1:end))])


		% check for convergence
		if (max(abs(y-yold))<opt.abstol)
			break;
		end
		iter = iter+1;
		if (iter > opt.maxiter)
			error('no convergence');
		end

	end % while 1

% TODO what is with the second variable? linear differences?
function apply_bc(id,bc)
	switch (bc.type)
	case {'level','surface'}
		Mzz(id,:)   = 0;
		Mzu(id,:)   = 0;
		Mzz(id,id)  = 1;
		rz(id)      = bc.val;
	case {'velocity'}
		Muz(id,:)   = 0;
		Muu(id,:)   = 0;
		Muu(id,id)  = 1;
		ru(id)      = bc.val;
	case {'discharge'}
		% Q = d*u
		Muz(id,:)  = 0;
		Muu(id,:)  = 0;
		Muu(id,id) = d(id);
		ru(id)     = bc.val;
	otherwise
		error('here');
	end % switch
end % apply_bc

function ds = ds(x)
	ds = diag(sparse(x));
end

end % function river_tide
