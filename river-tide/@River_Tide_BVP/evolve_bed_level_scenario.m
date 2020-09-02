% Mon 24 Aug 11:25:05 +08 2020
% Karl Kastner, Berlin
%
%% shortcut function for batch simulation runs
%
function [t,zb] = evolve_bed_level_scenario(obj ...
		 , z10    ...	% [1]     amplitude of incoming wave
		 , pz1r   ...	% [1]     reflected wave factor
		 , omega  ...	% [rad/s] anguluar frequency of tide
		 , Qmin   ...	% [m^3/s] discharge at inflow bc
		 , Qmax   ...	% [m^3/s] discharge at inflow bc
		 , iorder ...
		 , S0     ...	% [1]	  upstream slope
		 , d_mm   ...
		 , L      ...	% [m]     doamin length
		 , pL     ...	% [1]     relative location of bifurcation as seen
		          ...	%         from mouth with respect to domain length
		 , w0     ...	% [m]     width of upstream channel
		 , Cd     ...
		 , ruleQ  ...	%  -      rule for sediment division
		 , bif_ignore_rt    ...	%  -      rule for sediment division
		 , bif_stokes_order ... %
		 , dx     ...	% [m]     spatial discretization step
		 , Ti     ...		% [s]     simulated time span
		 , cfl    ...	% [1]     cfl condition
		 , scheme ... 	%  -      numerical scheme
		 , opt    ...
	)
	if (nargin()<18)
		opt              = struct();
	end

	Q0 = formative_discharge(Qmin,Qmax);

	obj.set_width(w0);
	obj.set_cd(Cd);

	h0 = normal_flow_depth(Q0,w0(1),Cd,S0,'cd');
	L0 = h0/S0;

	obj.omega = omega;

	% domain size
	Xi = [   L0*pL,L0;
	             0,L0*pL;
	             0,L0*pL];

%	meta = river_tide_test_metadata();
	
	opt.maxiter      = 100;
	opt.sopt.maxiter = 100;
	obj.opt.dischargeisvariable = true;
	obj.opt.imethod = 'spline';
	obj.opt = copy_fields(opt,obj.opt);

	bifurcation = Bifurcation();
	bifurcation.division_rule    = @bifurcation.sediment_division_geometric;
	bifurcation.opt.stokes_order  = bif_stokes_order; %obj.opt.stokes_order;
	bifurcation.opt.ignore_rt = bif_ignore_rt;
	obj.bifurcation = bifurcation;

	% perturbation of the bed level
	dz = [0,+0.1,-0.1]*h0;

	% initial bed level of channel
	zb        = @(cdx,x) -h0 + S0*x + dz(cdx);
	obj.set_zb(zb);

	obj.sediment.d_mm = d_mm;

	bc           = struct();
	bc_Qs        = struct();

	nx = [];
	for cdx=1:length(w0)
		nx(cdx) = round(diff(Xi(cdx,:))/dx);
	
		% downstream condition
		% mean sea level
		if (1==cdx)
			bc(1,1,cdx).var = '';
			bc(1,1,cdx).rhs = [];
		else
			bc(1,1,cdx).var = 'z';
			bc(1,1,cdx).rhs = 0;
			% Dirichlet condition
			bc(1,1,cdx).p   = 1;
		end
	
		% upstream condition
		% river discharge
		if (1==cdx)
			bc(2,1,cdx).var = 'Q';
			%bc(2,1,cdx).rhs = Q0;
			bc(2,1).Qseason = [Qmin,Qmax];
		else
			bc(2,1,cdx).var = '';
			bc(2,1,cdx).rhs = [];
		end
	
		% wave entering from left
		if (1==cdx)
			bc(1,2,cdx).var = '';
			bc(1,2,cdx).rhs = [];
		else
			bc(1,2,cdx).var = 'z';
			bc(1,2,cdx).rhs = z10;
			bc(1,2,cdx).p   = [1,0];
			bc(1,2,cdx).q   = [1,pz1r];
		end
	
		% wave entering from right
		if (1 == cdx)
			bc(2,2,cdx).var = 'z';
			bc(2,2,cdx).rhs =   0;
			bc(2,2,cdx).p   = [1,0];
			bc(2,2,cdx).q   = [0,1];
		else
			bc(2,2,cdx).var = '';
			bc(2,2,cdx).rhs = [];
		end

		if (1 == cdx)
			bc_Qs(1,cdx).rhs = [];
			bc_Qs(1,cdx).p   = 0;
			bc_Qs(2,cdx).rhsfun = @Qsfun; 
			bc_Qs(2,cdx).p   = 1;
		else
			bc_Qs(1,cdx).rhs = 0;
			bc_Qs(1,cdx).p   = 1;
			bc_Qs(2,cdx).rhs = [];
			bc_Qs(2,cdx).p   = 1;
		end
	end % for cdx (each channel)

	obj.bc    = bc;
	obj.bc_Qs = bc_Qs;

	function [cid,eid,p] = jfun()
		cid = [1,2,3];
		eid = [1,2,2];
		p   = [1./w0];
	end

	hydrosolver = BVPS_Characteristic();
	hydrosolver.xi = Xi;
	hydrosolver.nx = nx;

	morsolver        = Time_Stepper();
	morsolver.Ti     = [0,Ti];
	morsolver.cfl    = cfl;
	morsolver.scheme = scheme;

	obj.hydrosolver    = hydrosolver;
	obj.morsolver    = morsolver;

	obj.junction_condition = {@jfun};
	obj.junction_Qs = {@jfun};

%	obj.check_arguments();
	[t,zb] = obj.evolve_bed_level();

function Qs_ = Qsfun(t,Q0_)
	h0_ = normal_flow_depth(Q0_,w0(1),Cd,S0,'Cd');
	Qs_ = total_transport_engelund_hansen(drag2chezy(Cd),obj.sediment.d_mm,Q0_./(h0_*w0(1)),[],w0(1));
end
	
end % evolve_bed_level_scenario

