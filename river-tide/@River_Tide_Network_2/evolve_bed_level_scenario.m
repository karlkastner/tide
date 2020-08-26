% Mon 24 Aug 11:25:05 +08 2020
function [t,zb] = evolve_bed_level_scenario(obj...
		 , z10  ...	% [1]     amplitude of incoming wave
		 , pz1r ...	% [1]     reflected wave factor
		 , omega ...	% [rad/s] anguluar frequency of tide
		 , Q0 ...	% [m^3/s] discharge at inflow bc
		 , S0 ...	% [1]	  upstream slope
		 , d_mm ...
		 , L ...	% [m]     doamin length
		 , pL ...	% [1]     relative location of bifurcation as seen
		     ...	%         from mouth with respect to domain length
		 , w0 ...	% [m]     width of upstream channel
		 , Cd ...
		 , ruleQ ...	%  -      rule for sediment division
		 , ruleQs ...	%  -      rule for sediment division
		 , dx ...	% [m]     spatial discretization step
		 , Ti ...		% [s]     simulated time span
		 , cfl ...	% [1]     cfl condition
		 , scheme ... 	%  -      numerical scheme
		 , opt ...
	)
	if (nargin()<18)
		opt              = struct();
	end

% [kg/s]  sediment transport at inflow bc

	% river discharge
	%Q0   = -10;
	%z10  =   1;

%	h0  = 10;
%	w0  = [1,1/3,2/3];
%	w0 = [1,1/2,1/2];
	%w0  = [1,1/2,1/2];
%	nx  = 50;
%	Cd = 2.5e-3;
%	Cd_ = 2.5e-3;

	h0 = normal_flow_depth(Q0,w0(1),Cd,S0,'cd');

	L_ = h0/S0;
%	S0  = -normal_flow_slope(Q0,h0,w0(1),drag2chezy(Cd));
%	S0_  = -normal_flow_slope(Q0_,h0,w0(1),drag2chezy(Cd_));
%	L = h0/S0_;

	% domain size
	Xi = [      L_*pL,L;
	             0,L_*pL;
	             0,L_*pL];

%	meta = river_tide_test_metadata();
% meta.opt;
	
	opt.maxiter      =  100;
	opt.sopt.maxiter = 100;

	rt = River_Tide();
	obj.rt = rt;

	% perturbation of the bed level
	dz = [0,+0.1,-0.1]*h0;
	for idx=1:length(w0)
		nx = round(diff(Xi(idx,:))/dx);
	
		opt.nx = nx; % round(nx*px(idx));
	
		% width of channel
		wfun      = @(x) w0(idx)*ones(size(x));
	
		% drag/friction coefficient
		cdfun     = @(x)  Cd*ones(size(x));
	
		% bed level of channel
		zbfun     = @(x) -h0 + S0*x + dz(idx);
	
		% base frequency
		%T         = Constant.SECONDS_PER_DAY;
		%omega     = 2*pi/T;
	
		bc        = struct();
	
		% mean sea level
		if (1~=idx)
			bc(1,1).var = 'z';
			bc(1,1).rhs = 0;
	
			% Dirichlet condition
			bc(1,1).p   = 1;
		else
			bc(1,1).var = '';
			bc(1,1).rhs = [];
		end
	
		% river discharge
		if (1==idx)
			bc(2,1).var = 'Q';
			bc(2,1).rhs = Q0;
		else
			bc(2,1).var = '';
			bc(2,1).rhs = [];
		end
	
		% wave entering from left
		if (1~=idx)
			bc(1,2).var = 'z';
			bc(1,2).rhs = z10;
			bc(1,2).p   = [1,0];
			bc(1,2).q   = [1,pz1r];
		else
			bc(1,2).var = '';
			bc(1,2).rhs = [];
		end
	
		% wave entering from right
		if (1 == idx)
			bc(2,2).var = 'z';
			bc(2,2).rhs =   0;
			bc(2,2).p   = [1,0];
			bc(2,2).q   = [0,1];
		else
			bc(2,2).var = '';
			bc(2,2).rhs = [];
		end
	
		obj.rt(idx) = River_Tide( ...
					   'fun.zb',      zbfun ...
					 , 'fun.cd',      cdfun ...
					 , 'fun.width',   wfun ...
					 , 'omega',       omega ...
					 , 'opt',         opt ...
					 , 'Xi',          Xi(idx,:) ...
					);
		obj.rt(idx).bc = bc;
	
		bc_Qs        = struct();
		obj.rt(idx).sediment.d_mm = d_mm;
		if (1==idx)
			bc_Qs(1).p   = 0;
			bc_Qs(1).val = NaN;
			bc_Qs(2).p   = 1;
			Qs = total_transport_engelund_hansen(drag2chezy(Cd), ...
							     d_mm, ...
							     Q0./(h0*w0(idx)), ...
							     h0,w0(idx));
			bc_Qs(2).val = Qs; 
		else
			bc_Qs(1).p   = 0;
			bc_Qs(1).val = 0;
			bc_Qs(2).p   = 1;
			bc_Qs(2).val = NaN;
		end
		obj.rt(idx).bc_Qs    = bc_Qs;
	end % for idx (each branch)

	function [cid,eid,p] = jfun()
		cid = [1,2,3];
		eid = [1,2,2];
		p   = [1./w0];
	end

	%secpyear         = Physics.SECONDS_PER_YEAR;
	morsolver        = Time_Stepper();
	morsolver.Ti     = [0,Ti];
	morsolver.cfl    = cfl;
	morsolver.scheme = scheme;
%	morsolver.scheme = 'leapfrog-trapezoidal';

	switch (ruleQs)
	case {'geometric'}
		% nothing to do, default
	case {'geometric-nort'}
		obj.opt.ignorertfordivision = true;
	otherwise
		error('here');
	end
	
	obj.rt(1).morsolver    = morsolver;

	%rtn = River_Tide_Network_2(rt);
	obj.junction_condition = {@jfun}
	obj.junction_Qs = {@jfun}
	%rtn.opt.ignorertfordivision = true;

%	obj.rt = rt;
	
	[t,zb] = obj.evolve_bed_level();
end % evolve_bed_level_scenario

