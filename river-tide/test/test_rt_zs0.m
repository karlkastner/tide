% Sun 11 Mar 18:59:47 CET 2018

function test_rt_zs0(meta)
	addpath ../river-tide-analysis/
	if (nargin() < 1)
		meta = river_tide_metadata;
	end
	oname     = meta.filename.rt_experiment;
%	model_str = meta.rt_model_str;

	% estuary properties parameters

	range.S0 = [-6, -3];

	% number of parameter sets per decade
	n.S0     = 3;
	S0       = [0,logspace(range.S0(1),range.S0(2),(diff(range.S0))*n.S0+1)];
	S0       = 16/3e5;

	range.cd = log10(1e-3)*[1 1]; %-4, -2];
	n.cd     = 1;
	cd       = logspace(range.cd(1),range.cd(2),(diff(range.cd))*n.cd+1)
	cd = 2.5e-3;

	range.h0 = log10(16)*[1 1]; %[ 0.5, 1.5];
	n.h0     = 1;
	zb_downstream = -logspace(range.h0(1),range.h0(2),(diff(range.h0))*n.h0+1)
	h0 = 16;

	% boundary conditions
	range.z1 = log10(1)*[1 1];
	n.z1     = 1;
	%z1_downstream  = logspace(range.z1(1),range.z1(2),(diff(range.z1))*n.z1+1);
	z1_downstream  = 0;

	range.Q0 = log10(10)*[1,4];
	n.Q0     = 1;
	Q0       = [logspace(range.Q0(1),range.Q0(2),(diff(range.Q0))*n.Q0+1)];
	Q0 = [1e3,1e4];

	g             = Constant.g;
	W0            = 500;

	omega1        = 2*pi/86400;
	omega         = omega1;

	q = [1 0];

	% options for swe
	% opt.aa = 1;
	
	% computational parameter
	opt.precomputeh0 = false;

%	X              = [0 1e7];
	X              = [0 1e6];
%	opt.nx         = 1024;
%	opt.nx		= 10;
	% stretch
	opt.xs	       = 100;
	opt.nf         = 1;
	opt.sopt.relaxation = 0.5;
	opt.sopt.maxiter    = 100;
	opt.model_str  = 'wave';
	opt.solver     = @bvp2fdm;
	opt.nx         = 1024;
	%opt.solver     = @bvp2c;

	val_C = {      {X} ...
		       , Q0 ...
		       , W0 ...
		       , S0 ...
		       , z1_downstream ...
		       , cd ...
		       , zb_downstream ...
		       , omega ...
		       , q ...
		       , opt }

	rtmap = RT_map(oname);
	rtmap.init();
	rtmap.recompute = true;
	% TODO pass arguments
	rtmap.run(val_C,true)
	rtmap.save()

end % test_rt_zs0

