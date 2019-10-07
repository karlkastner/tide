% Wed  2 May 19:22:47 CEST 2018
% TODO, exclude reflected wave
	addpath ../river-tide-analysis/
	meta 	  = river_tide_metadata;
	iname     = meta.filename.rt_experiment;
	model_str = meta.rt_model_str;
	if (~exist('pflag','var'))
		pflag = 0;
	end
	fflag = pflag;
	ps = 2.5;

	if (~exist('rtmap','var')) % || reload)
		rtmap = RT_map(iname);
		%load(iname,'rt');
		rtmap.init();
	end

	p = 1/4;
	g = Constant.gravity;
	z1_downstream = sqrt(eps);
	omega   = 2*pi/86400;
	cd 	= 2.5e-3;
	h0	= 15;
	L       = 3e5;
	dzb_dx  = h0/L;
	dzb_dx_ = dzb_dx;
	W       = 500;
	C       = sqrt(g/cd);
	Qn      = normal_flow_discharge(h0,W,C,dzb_dx)

	Q0      = Qn;
	X       = [0,10*L];
	q       = [0,1];
	opt     = struct();

	opt.nx  = 1024;
	opt.xs	= 16;

	opt.model_str = 'wave';
	opt.solver     = @bvp2fdm;
	opt.o2 = true;
%	opt.o2 = false;

	% z1 admittance pre-legend
	splitfigure([2 2],[2 1],fflag);
	cla();
	subplot(2,2,1)
	plot(NaN(length(Q0)));
	hold on
	set(gca,'colororderindex',1);
	cmap = colormap('lines');



	% constant cross section along channel

	% TODO, make rtmap arguments a struct
	rt = rtmap.fun( ...
					{X}, ...
					Q0, ...
					W, ...
					dzb_dx, ...
					z1_downstream, ...
					cd, ...
					-h0, ...
					omega, ...
					{q}, ...
					opt ...
			);
	rt.opt.hmode = 'predefined';
	% TODO introduce extract_unknowns and stack_unknowns
	rt.opt.o2    = false;

	for idx=1:2

	if (1 == idx)
		rt.fun.z0 = @(x) zeros(size(x));
		rt.fun.zb = @(x) -h0*ones(size(x));
	else
		% variable cross section along channel
		Lh = 1e5/4;
		rt.fun.z0 = @(x) zeros(size(x));
		rt.fun.zb = @(x) -h0*exp(-x/Lh);	
	end
	rt.init([0,z1_downstream], Q0, X);
	rt.solve();

	x = rt.x;
	az1=rt.admittance('z',1);
	aq1=rt.admittance('q',1);
	kz1 = rt.wave_number('z',1);
	kq1 = rt.wave_number('q',1);
	[k10, kz1a, kq1a, dk, obj] = rt.wave_number_approximation();
	id = 1;
	dk = dk(:,id);
	kq1a = kq1a(:,id);
	kz1a = kz1a(:,id);
	k10 = k10(:,id);

	fdx = find(az1<1e-2,1,'first');
	x0 = x(fdx);

	figure(idx);
	clf();
	subplot(2,2,1)
	plot(x,[imag(kz1),real(kz1)],'linewidth',1);
	hold on;
	set(gca,'colororderindex',1)
	plot(x,[imag(kz1a), real(kz1a)],'--','linewidth',2)
	set(gca,'colororderindex',1)
	plot(x,[imag(k10), real(k10)],'-.','linewidth',1)
	%xlim([0,x(end)/2])
	legend('im','re','ima','re');
	xlim([0,x0])
	ylabel('k_z');

	subplot(2,2,2);
	plot(x,[imag(kq1),real(kq1)],'linewidth',1);
	hold on;
	set(gca,'colororderindex',1)
	plot(x,[imag(kq1a), real(kq1a)],'--','linewidth',2)
	set(gca,'colororderindex',1)
	plot(x,[imag(k10), real(k10)],'-.','linewidth',1)
	%xlim([0,x(end)/2])
	legend('im','re','im appr','re appr','im k10','re k10');
	xlim([0,x0])
	ylabel('k_q');

	subplot(2,2,3);
	plot(x,rt.h0,'linewidth',1);
	xlim([0,x0])
	ylabel('h0');

	subplot(2,2,4);
	plot(x,[az1 aq1],'linewidth',1);
	xlim([0,x0])
	ylabel('|z_1|');

	figure(3);
	subplot(2,2,1)
	plot(x,[imag(k10) -imag(dk)]); %, real(dk)./real(k10)])
	xlim([0,x0])
	subplot(2,2,2)
	plot(x,[imag(dk)./imag(k10), real(dk)./real(k10)])
	xlim([0,x0])
	

	
	[norm(select(kz1-kz1a,1:fdx)), 	norm(select(kz1-k10,1:fdx))]
	[norm(select(kq1-kq1a,1:fdx)), 	norm(select(kq1-k10,1:fdx))]
	
	
	end
