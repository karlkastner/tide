% Mon  2 Oct 12:27:30 CEST 2017
%% compute river tide for a scenario with specific boundary conditions and store it in the hash,
%% or retrive the scenario, if it was already computed
function [out, key, obj] = fun(obj ...
		, Xi ... 	% domain extend (river start and end coordiante)
		, width ... 	% channel width as scalar or function of x
		, cd ...	% drag coefficient as scalar or function x
		, zb ...	% bed level as scalar or function of x
		, omega ...	% angalar frequency of tide
		, bc0l_var, bc0l_val, bc0l_p ...	% mean flow bc at mouth
		, bc0r_var, bc0r_val, bc0r_p ...	% mean flow bc upstream
		, bc1l_var, bc1l_val, bc1l_p, bc1l_q ... % main species bc at mouth
		, bc1r_var, bc1r_val, bc1r_p, bc1r_q ... % main species bc upstream
		, bc2l_var, bc2l_val, bc2l_p, bc2l_q ... % overtide bc at mouth
		, bc2r_var, bc2r_val, bc2r_p, bc2r_q ... % overtide bc upstream
		, opt ... % additional options
		)

	if (nargin()<10)
		opt = [];
	end

	if (isempty(obj.cmap))
		obj.init();
	end

	g     = Constant.g;
	Xi    = Xi{1};
	L     = Xi(2)-Xi(1);
	width     = width{1};
	cd    = cd{1};
	zb    = zb{1};
	
	if (~isempty(bc2l_var))
		bc2l = {bc2l_var{1}, ...
			bc2l_val{1}, ...
			bc2l_p{1}(1), ...
			bc2l_p{1}(2), ...
			bc2l_q{1}(1), ...
			bc2l_q{1}(2)};
	else
		bc2l = repmat({},6,1);
	end
	if (~isempty(bc2r_var))
		bc2r = {bc2r_var{1}, ...
			bc2r_val{1}, ...
			bc2r_p{1}(1), ...
			bc2r_p{1}(2), ...
			bc2r_q{1}(1), ...
			bc2r_q{1}(2)};
	else
		bc2r = repmat({},6,1);
	end

	if (isa(zb,'function_handle'))
		zb_str = func2str(zb);
	else
		zb_str = num2str(zb);
	end
	if (isa(width,'function_handle'))
		w_str = func2str(width);
	else
		w_str = num2str(width);
	end
	if (isa(cd,'function_handle'))
		cd_str = func2str(cd);
	else
		cd_str = num2str(cd);
	end

	key = obj.key(  ...	
			...  % opt.model_str, ...
			...  % func2str(opt.solver), ...
			zb_str, ...
			w_str, ...
			cd_str, ...
			omega, ...		
			bc0l_var{1}, ...		% left end
			bc0l_val{1}, ...
			bc0l_p{1}(1), ...
			bc1l_var{1}, ...
			bc1l_val{1}, ...
			bc1l_p{1}(1), ...
			bc1l_p{1}(2), ...
			bc1l_q{1}(1), ...
			bc1l_q{1}(2), ...
			bc2l{:}, ...
			bc0r_var{1}, ...		% right end
			bc0r_val{1}, ...
			bc0r_p{1}(1), ...
			bc1r_var{1}, ...
			bc1r_val{1}, ...
			bc1r_p{1}(1), ...
			bc1r_p{1}(2), ...
			bc1r_q{1}(1), ...
			bc1r_q{1}(2), ...
			bc2r{:}, ...
			L,	 ...
			opt.nx,	 ...
			opt.xs ...
		);

	% test if simulation was not yet run for the parameter pair
	if (isKey(obj.cmap,key) && ~obj.recompute)
		out = obj.cmap(key);
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
			out = rt_quasi_stationary_trigonometric(Xi,W0,Q0_,cd,omega,zbfun,bc,opt);
		case {'wave','swe'}
			% x         = linspace(Xi(1),Xi(2))';

			hydrosolver    = BVPS_Characteristic();
			hydrosolver.xi = Xi;
			hydrosolver.nx = opt.nx;
			hydrosolver.opt.xs = opt.xs;

			opt.dischargeisvariable = true;

			out = River_Tide_BVP( ...
				   'zb',        zb ...
				 , 'cd',        cd ...
				 , 'width',     width ...
				 , 'omega',         omega ...
				 , 'opt',           opt ...
				 , 'hydrosolver',   hydrosolver ...
				);

			bc            = out.bc;
		
			% boundary condition mean component (left end)
			bc(1,1).var   = bc0l_var{1};
			bc(1,1).rhs   = bc0l_val{1};
			bc(1,1).p     = bc0l_p{1};
			% right end
			bc(2,1).var   = bc0r_var{1};
			bc(2,1).rhs   = bc0r_val{1};
			bc(2,1).p     = bc0r_p{1};

			% boundary condition main tidal component
			bc(1,2).var   = bc1l_var{1};
			bc(1,2).rhs   = bc1l_val{1};
			bc(1,2).p     = bc1l_p{1}; 
			bc(1,2).q     = bc1l_q{1};
			bc(2,2).var   = bc1r_var{1};
			bc(2,2).rhs   = bc1r_val{1};
			bc(2,2).p     = bc1r_p{1}; 
			bc(2,2).q     = bc1r_q{1};

			if (~isempty(bc2l_var))
			% bc of even overtide
				bc(1,3).var   = bc2l_var{1};
				bc(1,3).rhs   = bc2l_val{1};
				bc(1,3).p     = bc2l_p{1}; 
				bc(1,3).q     = bc2l_q{1}; 
			end
			if (~isempty(bc2r_var))
				bc(2,3).var   = bc2r_var{1};
				bc(2,3).rhs   = bc2r_val{1};
				bc(2,3).p     = bc2r_p{1}; 
				bc(2,3).q     = bc2r_q{1}; 
			end
		
			out.bc = bc;

			out.init();

			out.solve();

			obj.cmap(key) = out;
		otherwise
			error('unimplemented solver');
		end
	end % if ~iskey
end % River_Tide_Map/fun

