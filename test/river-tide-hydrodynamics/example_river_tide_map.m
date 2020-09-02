% Thu 10 Oct 12:52:41 PST 2019
% exmaple river tide

	tid   = 1000;
	pflag = 1;
	name  = 'Example River\_Tide\_Map';

	% river discharge
	Qu          = -10;
	Q0          = Qu;
	Q0          = linspace(0.1,2,5)*Qu;
	% width of channel
	w0          = 1;
	wfun        = @(x)   w0*ones(size(x));
	% drag/friction coefficient
	cD          = 2.5e-3;
	cdfun       = @(x)  cD*ones(size(x));
	% bed level of channel
	h0          = 10;
	S0          = normal_flow_slope(-Qu,h0,w0,drag2chezy(cD));
	zbfun       = @(x) -h0 + S0*x; %*ones(size(x));
	% mean component
	bc          = struct();
	% mean sea level
	bc(1,1).var = 'z';
	bc(1,1).rhs = 0;
	% Dirichlet condition
	bc(1,1).p   = 1;
	% river discharge
	bc(2,1).var = 'Q';
	bc(2,1).rhs = -Q0;
	% Dirichlet condition
	bc(2,1).p   = 1;

	% wave entering from left
	bc(1,2).var = 'z';
	z10         = sqrt(eps);
	bc(1,2).rhs = z10;
	bc(1,2).p   = [1,0];
	bc(1,2).q   = [1,0];

	% wave entering from right / reflected wave
	bc(2,2).var = 'z';
	bc(2,2).rhs =   0;
	bc(2,2).p   = [1,0];
	bc(2,2).q   = [0,1];

	% base frequency
	T           = Constant.SECONDS_PER_DAY;
	omega       = 2*pi/T;
	% domain size
	Xi          = [0,1e6];
	% model for river tide
	opt.model_str = 'wave';
	% solver of boundary value problem
	opt.solver  = @bvp2c;
	% number of points along channel
	opt.nx      = 200;
	% change of distance between points along channel 
	opt.xs      = 1; 

	% solve with model
	rt_map = River_Tide_Map('./example.mat');
	rt_map.run({{Xi}  ... % -Q0, 
			, {wfun}, {cdfun}, {zbfun}, omega ...
			, {bc(1,1).var}, {bc(1,1).rhs}, {bc(1,1).p} ...
			, {bc(2,1).var}, num2cell(bc(2,1).rhs), {bc(2,1).p} ...
			, {bc(1,2).var}, {bc(1,2).rhs}, {bc(1,2).p}, {bc(1,2).q} ...
			, {bc(2,2).var}, {bc(2,2).rhs}, {bc(2,2).p}, {bc(2,2).q} ...
			, opt});
	
	% fetch and plot
	z0 = [];
	z1 = [];
	leg_C = {};
	leg_C = arrayfun(@(x) sprintf('q_0 = %5.2f m^2/s',x), Q0,'uniformoutput',false);
	for idx=1:length(Q0)
	[rt]  = rt_map.fun({Xi}  ... % -Q0, 
			, {wfun}, {cdfun}, {zbfun}, omega ...
			, {bc(1,1).var}, {bc(1,1).rhs}, {bc(1,1).p} ...
			, {bc(2,1).var}, {-Q0(idx)}, {bc(2,1).p} ...
			, {bc(1,2).var}, {bc(1,2).rhs}, {bc(1,2).p}, {bc(1,2).q} ...
			, {bc(2,2).var}, {bc(2,2).rhs}, {bc(2,2).p}, {bc(2,2).q} ...
			, opt);
	z0(:,idx) = rt.z(0); 
	z1(:,idx) = rt.z(1)); 
	end % for idx

	namedfigure(1,name);
	clf();
	subplot(2,2,1);
	plot(rt.x,z0);
	legend('location','northwest',leg_C{:});
	xlabel('x / m');
	ylabel('z_0 / m');
	
	subplot(2,2,2);
	plot(rt.x,abs(z1));
	xlabel('x / m');
	ylabel('|z_1| / m');
	legend(leg_C{:});

	% for later reuse of computation call:
	% rt_map.save();

