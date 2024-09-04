% Fri  2 Sep 15:12:52 CEST 2022
function rt = hydrodynamic_scenario_from_table(rt_map,file_str,id,opt)

	tab = readtable(file_str);

	fdx   = find(tab.id == id);
	
	% note : variables are expanded here, to allow use in eval calls

	description  = tab(fdx,:).description{1};

	% initial surface elevation
	%z00 = 0;

	% surface elevation at channel mouth
	z10 = tab.z10(fdx);

	z20 = tab.z20(fdx);

	% tidal surface elevation
	zs  = [0,z10,z20]

	% river discharge
	Q0 = tab.Q0(fdx);
	if (iscell(Q0))
		Q0 = eval(Q0{1});
	end
	if (~isscalar(Q0))
		Qseason = Q0;
		Q0 = formative_discharge(Qseason(1),Qseason(2),'chezy');
	else
		Qseason = [];
	end
	%w00 = tab.w00(fdx);

	% width of channel
	w0  = eval(tab(fdx,:).w0{1});

	% drag/friction coefficient
	Cd = tab.Cd(fdx);

	% depth of channel at mouth
	h0 = tab.h0(fdx); 

	% width of channel at river mouth
	w00 = w0(0);

	% slope of channel bed at river mouth
	S0  = -normal_flow_slope(Q0,h0,w00,drag2chezy(Cd));

	% bed level of channel
	zb        = eval(tab(fdx,:).zb{1});

	% funademental period in days
	T1_d      = tab.T1(fdx);
	% fundamental period in seconds
	T1        = T1_d*Constant.SECONDS_PER_DAY;
	% fundamental angular frequency
	omega     = 2*pi/T1;

	T2_d      = tab.T2(fdx);
	if (isfinite(T2_d))
		omega = [0, omega(1), 2*pi/(Constant.SECONDS_PER_DAY*T2_d)];
	end

	% seasonal / low frequent river period
	if (any(ismember(tab.Properties.VariableNames,'Ts')))
		Ts = tab.Ts(fdx);
		if (isnan(Ts))
			Ts = [];
		end
	else
		Ts = [];
	end

	% length of computational domain
	Lx = tab.Lx(fdx);

	% reflection coefficient at end of boundary
	ql = tab.ql(fdx);
	qr = tab.qr(fdx);

	nx = tab.nx(fdx);

	rt = hydrodynamic_scenario(rt_map,zs,ql,qr,zb,Q0,w0,Cd,omega,Lx,opt,Ts,Qseason,nx,description);

end


