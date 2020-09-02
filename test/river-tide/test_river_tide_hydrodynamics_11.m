% Wed  9 Oct 15:23:10 PST 2019
function [out,rt] = river_tide_test_11(rt_map,pflag)
	if (nargin()<2)
		pflag = 1;
	end
	out.id  = 11;
	out.name = 'D1 and D2 comparison';
	
	% surface elevation at mouth
	z00 = 0;
	z20_ = 0.1;
	z10 = [0,z20_];
	z20 = [z20_,0];	
	
	% base frequency
	% period (in days)
	T   = [1,0.5];
	T = T*Constant.SECONDS_PER_DAY;
	omega     = 2*pi./T;

	% river discharge
	Q0        =  -10;
	% width of channel
	w0        = 1;
	% drag/friction coefficient
	Cd        = 2.5e-3;
	% depth of channel
	h0        = 10;
	% slope of channel bed
	S0        = -normal_flow_slope(Q0,h0,w0,drag2chezy(Cd));
	% bed level of channel
	zb     = @(x) -h0 + S0*x;
	% domain length
	Lx = 1e6;

	meta = river_tide_test_metadata();
	opt = meta.opt;
	opt.oflag   = [true(1,3)];

	% reflection coefficient at right end of boundary
	qr = 0;
	for idx=1:2
		rt(idx) = hydrodynamic_scenario(rt_map,[0,z10(idx),z20(idx)],qr,zb,Q0,w0,Cd,omega(idx),Lx,opt);
	end % for idx

	z = ([rt(1).z(2),rt(2).z(1)]);
	res = z(:,1)-z(:,2);
	rmse = max(abs(res));
	result = (rmse > 0.01*z20_);

	% dummy value
	rmse(2) = NaN;
	out.rmse   = rmse;
	out.result = result;


	if (pflag)
		figure(100+tid);
		clf();
		subplot(2,2,1)
		plot(rt(1).x,abs(z(:,1)));
		hold on
		plot(rt(1).x,abs(z(:,2)),'--');
		legend('T_{base} = 1','T_{base} = 1/2');
		ylabel('|z_2|');
	
		subplot(2,2,2)
		z = ([rt(1).z(2),rt(2).z(1)]);
		plot(rt(1).x,angle(z(:,1)));
		hold on
		plot(rt(1).x,angle(z(:,2)),'--');
		ylabel('arg(z_2)');
	
		subplot(2,2,3)
		plot(rt(1).x,abs([rt(1).out.z(:,2:3)]));
	end % if pflag

	% compare to analytical solution
%	g = Constant.gravity;
%	c = sqrt(g*h0);
%	k = omega/c;
%	x = rt.x;

%	bw = Backwater1D();
%	nn = opt.nx;
%	[x_, h_, z0_] = bw.solve(-Q0,0,drag2chezy(Cd),wfun,zbfun,0,Xi);
%	[x_, h_, z0_] = bw.solve_analytic(-Q0,drag2chezy(Cd),w0,S0,h0,nn);
%	z0  = interp1(x_,z0_,rt.x,'linear','extrap');
%	z0 = S0*x;
%	z0t = rt.z(0) - z0;
%	z2 = rt.z_(:,3);
%
%	z2_ = rt.even_overtide_analytic(z10);
%	%z0_ = interp1(x_,h_,rt.x,'spline')+0*zbfun(rt.x);
%	% r = (1+1i)*sqrt(-Cd.*omega.*Q0/w0./(g*h0.^3));
%	% z = z10*exp(-r*x);
%
%	rmse(2)  = norm(z2_-z2);
%
%	% err ~ C*df^2/dx^2*dx^2, where C sufficiently large constant
%	nres_ = rms(cdiff(z2_,2))
%	fail = (rmse(1) > 0.05/z10) || (rmse(2) > 10*nres_);
%
%	zz = NaN(size(rt.out.z));
%	zz(:,3) = z2_;
%	river_tide_test_plot(out.id,rt,zz,out.name,pflag);
%	
end % river_tide_test_09

