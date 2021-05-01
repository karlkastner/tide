% Wed  9 Oct 15:23:10 PST 2019
function [out, rt, d3d] = test_river_tide_hydrodynamics_11(rt_map,pflag)
	meta = test_river_tide_metadata();
	if (nargin()<1)
		rt_map      = River_Tide_Hydrodynamics_Map(meta.mapname_str);
	end
	if (nargin()<2)
		pflag = 1;
	end
	tab      = readtable('test-river-tide.csv');
	out.id   = 11;
	fdx      = find(tab.id == out.id)
	out.name = tab(fdx,:).name{1};

	% surface elevation at mouth
	z00 = 0;
	z20_ = tab.z10(fdx);
	z10 = [0,z20_];
	z20 = [z20_,0];	
	
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

	% slope of channel bed
	S0        = -normal_flow_slope(Q0,h0,w00,drag2chezy(Cd));

	% bed level of channel
	zb        = eval(tab(fdx,:).zb{1});

	% base frequency
	% period (in days)
	T   = [1,0.5]*Constant.SECONDS_PER_DAY;
	omega     = 2*pi./T;

	% length of computational domain
	Lx = tab.Lx(fdx);

	% reflection coefficient at right end of boundary
	ql = tab.ql(fdx);
	qr = tab.qr(fdx);

	opt = meta.opt;
	opt.oflag   = [true(1,3)];

	for idx=1:2
		rt(idx) = hydrodynamic_scenario(rt_map,[0,z10(idx),z20(idx)],ql,qr,zb,Q0,w0,Cd,omega(idx),Lx,opt);
	end % for idx

	% generate d3d equivalent model for comparison
	d3dopt                = struct();
	d3dopt.Lc            = tab.Lc(fdx);
	d3dopt.bndisharmonic = true;
	folder = [meta.folder.d3d,num2str(out.id)];
	rt(1).generate_delft3d(folder,meta.param,meta.param_silent,d3dopt);
	[out.rmse_d3d, d3d] = test_rt_d3d_evaluate(rt(1),out.id,pflag);

	z = ([rt(1).channel(1).waterlevel(2),rt(2).channel(1).waterlevel(1)]);
	res = z(:,1)-z(:,2);
	rmse = max(abs(res));
	result = (rmse > 0.01*z20_);

	% dummy value
	rmse(2) = NaN;
	out.rmse   = rmse;
	out.result = result;


	if (pflag)
		figure(100+out.id);
		clf();
		subplot(2,2,1)
		plot(rt(1).channel(1),abs(z(:,1)));
		hold on
		plot(rt(1).channel(1),abs(z(:,2)),'--');
		legend('T_{base} = 1','T_{base} = 1/2');
		ylabel('|z_2|');
	
		subplot(2,2,2)
		z = ([rt(1).waterlevel(2),rt(2).waterlevel(1)]);
		plot(rt(1).channel(1),angle(z(:,1)));
		hold on
		plot(rt(1).channel(1),angle(z(:,2)),'--');
		ylabel('arg(z_2)');
	
		subplot(2,2,3)
		plot(rt(1).channel(1),abs([rt(1).channel(1).z(:,2:3)]));
	end % if pflag


end % river_tide_test_11

