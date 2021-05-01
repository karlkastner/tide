% Wed  9 Oct 15:23:10 PST 2019
function [out, rt, d3d] = test_river_tide_hydrodynamics_12(rt_map,pflag)
	meta = test_river_tide_metadata();
	if (nargin()<1)
		rt_map      = River_Tide_Hydrodynamics_Map(meta.mapname_str);
	end
	if (nargin()<2)
		pflag = 1;
	end
	reltol = 1e-2;

	tab = readtable('test-river-tide.csv');
	out.id  = 12;
	fdx = find(tab.id == out.id)
	out.name = tab(fdx,:).name{1};

	% surface elevation at mouth
	z10_ = tab.z10(fdx);
	z10 = [z10_,0];
	z1L = [0,z10_];

	% river discharge
	Q0 = tab.Q0(fdx);

	% width of channel at river mouth
	w00 = tab.w00(fdx);

	% width of channel
	w0  = eval(tab(fdx,:).w0{1});

	% drag/friction coefficient
	Cd = tab.Cd(fdx);

	% depth of channel
	h0 = tab.h0(fdx);

	for idx=1:2
	if (2 == idx)
		Q0 = -Q0;
	end

	% slope of channel
	S0         = -normal_flow_slope(Q0,h0,w00,drag2chezy(Cd))

	% length of computational domain
	Lx = tab.Lx(fdx);

	% bed level of channel
	if (1==idx)
		zb  = @(x) -h0 + S0*x;
	else
		
		zb  = @(x) -h0 + S0*(x-Lx);
	end
	bc        = struct();

	% mean sea level
	if (1 == idx)
	bc(1,1).var = 'z';
	bc(1,1).rhs = 0;
	% river discharge
	bc(2,1).var = 'Q';
	bc(2,1).rhs = Q0;
	else
	bc(1,1).var = 'Q';
	bc(1,1).rhs = Q0;
	% river discharge
	bc(2,1).var = 'z';
	bc(2,1).rhs = 0;
	end
	% Dirichlet condition
	bc(1,1).p   = 1;
	% Dirichlet condition
	bc(2,1).p   = 1;
	q = 0;

	% wave entering from left
	bc(1,2).var = 'z';
	bc(1,2).rhs = z10(idx);
	bc(1,2).p   = [1,0];
	bc(1,2).q   = [1,q];

	bc(1,3).var = 'z';
	bc(1,3).rhs = 0;
	bc(1,3).p   = [1,0];
	bc(1,3).q   = [1,0];

	% wave entering from right / reflected wave
	bc(2,2).var = 'z';
	bc(2,2).rhs = z1L(idx);
	bc(2,2).p   = [1,0];
	bc(2,2).q   = [q,1];

	bc(2,3).var = 'z';
	bc(2,3).rhs = 0;
	bc(2,3).p   = [1,0];
	bc(2,3).q   = [0,1];

	% base frequency
	% period (in days)
	T_d       = tab.T(fdx);
	T   = T_d*Constant.SECONDS_PER_DAY;
	omega     = 2*pi/T;

	Xi        = [0,Lx];

	% reflection coefficient at right end of boundary
	ql = tab.ql(fdx);
	qr = tab.qr(fdx);

	opt = meta.opt;
	opt.oflag   = [true(1,1)];

	% solve with model
	[rt(idx)]  = rt_map.fun({Xi} ... % Q0,
			, {w0}, {Cd}, {zb}, omega ...
			, {bc(1,1).var}, {bc(1,1).rhs}, {bc(1,1).p} ...
			, {bc(2,1).var}, {bc(2,1).rhs}, {bc(2,1).p} ...
			, {bc(1,2).var}, {bc(1,2).rhs}, {bc(1,2).p}, {bc(1,2).q} ...
			, {bc(2,2).var}, {bc(2,2).rhs}, {bc(2,2).p}, {bc(2,2).q} ...
			, {bc(1,3).var}, {bc(1,3).rhs}, {bc(1,3).p}, {bc(1,3).q} ...
			, {bc(2,3).var}, {bc(2,3).rhs}, {bc(2,3).p}, {bc(2,3).q} ...
			, opt);

	end

	% generate d3d equivalent model for comparison
	d3dopt                = struct();
	d3dopt.Lc            = tab.Lc(fdx);
	d3dopt.bndisharmonic = true;
	folder = [meta.folder.d3d,num2str(out.id)];
	rt(1).generate_delft3d(folder,meta.param,meta.param_silent,d3dopt);
	[out.rmse_d3d, d3d] = test_rt_d3d_evaluate(rt(1),out.id,pflag);


	zb = ([rt(1).channel(1).zb,flipud(rt(2).channel(1).zb)]);

	z0 = ([rt(1).channel(1).waterlevel(0),flipud(rt(2).channel(1).waterlevel(0))]);

	z1 = ([rt(1).channel(1).waterlevel(1),flipud(rt(2).channel(1).waterlevel(1))]);
	res = z1(:,1)-z1(:,2);
	rmse = rms(res);
	result = (rmse(1) > reltol*rms(z1(:,1)));

	% dummy value
	rmse(2) = NaN;
	out.rmse   = rmse;
	out.result = result;

	if (pflag)
		figure(1);
		clf();
		subplot(2,2,1)
		plot(rt(1).channel(1).x,[z0(:,1) zb(:,1)]);
		hold on
		plot(rt(1).channel(1).x,[z0(:,2) zb(:,2)],'--');
		ylabel('z_0');
	
		subplot(2,2,3)
		plot(rt(1).channel(1).x,abs(z1(:,1)));
		hold on
		plot(rt(1).channel(1).x,abs(z1(:,2)),'--');
		ylabel('|z_1|');
	
		subplot(2,2,4)
		plot(rt(1).channel(1).x,angle(z1(:,1)));
		hold on
		plot(rt(1).channel(1).x,angle(z1(:,2)),'--');
		ylabel('arg(z_1)');
	
		subplot(2,2,2)
		plot(rt(1).channel(1).x,[rt(1).channel(1).zb,rt(2).channel(1).zb],'k');
		hold on
		plot(rt(1).channel(1).x,[rt(1).channel(1).waterlevel(0),rt(2).channel(1).waterlevel(0)],'r--');
	end


end % river_tide_test_12

