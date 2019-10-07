% Mon  2 Oct 12:27:30 CEST 2017
%% compute river tide for a scenario with specific boundary conditions and store it in the hash,
%% or retrive the scenario, if it was already computed
function [out, key, obj] = fun(obj,Xi,Q0,W0,S0,z1_downstream,cd,zb_downstream,omega,q,opt)
	if (nargin()<10)
		opt = [];
	end
%	if (~isfield(opt,'model_str'))
%		opt.model_str = 'wave';
%	end
	g   = Constant.g;
	C   = sqrt(g/cd);
	Xi  = Xi{1};
	L   = Xi(2)-Xi(1);
	q   = q{1} 
	key = obj.key(...	
			opt.model_str, ...
			func2str(opt.solver), ...
			zb_downstream, ...
			S0, ...
			W0, ...
			cd, ...
			Q0, ...
			z1_downstream, ...
			omega,   ...
			q(1), ...
			q(2), ...
			L,	 ...
			opt.nx,	 ...
			opt.xs ...
		);

	% test if simulation was not yet run for the parameter pair
	if (isKey(obj.rt,key) && ~obj.recompute)
		out = obj.rt(key);
	else
		disp(['recomputing ', key]);
		obj.recflag = true;
%
		switch (opt.model_str)
		case {'odeset'}
			bc = struct();
			if (0)
				fzb = @(x) zb_downstream+(Xi(end)-x)*S0;
				bc.z0  = [0 0];
				%fzb(Xi(1))+normal_flow_depth 
				bc.cz1 = [0 z1_downstream]; % [0 1]
				bc.sz1 = [0 0.0];
			else
				fzb = @(x) zb_downstream + x*S0;
				Q0_     = -Q0;
				bc.z0  = [0 0];
				bc.cz1 = [z1_downstream 0]; % [1 0]
				bc.sz1 = [0 0];
			end
			out = rt_quasi_stationary_trigonometric(Xi,W0,Q0_,cd,omega,fzb,bc,opt);
		case {'wave','swe'}
			% flag     = [];
			bc     = struct();
			bc.p   = [0, 1];
			bc.q   = q;
			bc.rhs = -1i*W0*omega*z1_downstream;

			zbfun_ = @(x) zbfun(x,zb_downstream,S0);

			x         = linspace(Xi(1),Xi(2))';
			[zb_, s_] = zbfun_(x);

			% TODO make RT intelligent to accept scalars instead of functions
			out = River_Tide( ...
				   'fun.zb',  zbfun_ ...
				 , 'cdfun',  @(x) cd*ones(size(x)) ...
				 , 'wfun',   @(x) W0*ones(size(x)) ...
				 , 'omega',  omega ...
				 ... % , 'bc',     bc ...
				 , 'opt',    opt ...
				 ... %, 'flag', flag ...
				);
			out.bc(1,2).p   = bc.p;
			out.bc(1,2).q   = bc.q;
			out.bc(1,2).rhs = bc.rhs;

			out.init([0,z1_downstream], Q0, Xi);

			out.solve();

			obj.rt(key) = out;
		otherwise
			error('unimplemented solver');
		end
	end % if ~iskey

	function [zb, dzb_dx] = zbfun(x,zb_downstream,S0)
		p   = 0.5;
		L   = Xi(2)-Xi(1);
		zL  = 1;
		zL_ = zb_downstream + S0*L;

		zb     = zb_downstream + S0*x;
		dzb_dx = S0*ones(size(x));

		if (zL_ < zL)
			fdx = (x>p*L);
			n   = sum(fdx);
			if (n>0)
				% there is reflection wherever the bed slope changes not gradually,
				% so the bed level is made c1 continuous

				% scale x by L
			
				Lhs = [vander_1d(p,2)
				       vanderd_1d(p,1,2)
				       vander_1d(1,2)];
				rhs = [zb_downstream + p*S0*L;S0;zL];
				c   = Lhs \ rhs;
				Lhs = vander_1d(x(fdx)/L,2);
				zb(fdx) = Lhs*c;

				Lhs = vander_1d(x(fdx)/L,1);
				c   = [c(2);2*c(3)];
				dzb_dx(fdx) = Lhs*c;
			end
		end
	end % sbfun

end % River_Tide_Map/fun

