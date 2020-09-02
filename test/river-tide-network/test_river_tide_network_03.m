% Tue 18 Aug 09:27:08 +08 2020
function [out] = river_tide_network_test_02(rt_map,pflag)
	if (nargin()<2)
		pflag = 1;
	end
	out.id  = 2;
	out.name = 'two channels, reflection at closed end';

	% tidal amplitude
	z10 = 1;

	% river discharge
	Q0  = 0;

	% width of channels
	w0  = 1*[1,1/3];

	% drag/friction coefficient
	Cd = 0;

	% depth of channels
	h0  = 10;

	% slope of channel bed
	S0        = -normal_flow_slope(Q0,h0,w0(1),drag2chezy(Cd));

	% bed level of channel
	zb     = @(x) -h0 + S0*x;

	% base frequency
	T         = Constant.SECONDS_PER_DAY;
	omega     = 2*pi/T;

	% domain size
	L = 8e5;
	Xi        = L*[   0,0.5;
                        0.5,1];

	% number of grid points
	nx  = [50, 51];

	meta = river_tide_test_metadata();
	opt = meta.opt;
	opt.maxiter =  100;
	opt.sopt.maxiter = 100;

	bc        = struct();
	for idx=1:2

	% downstream boundary condition
	if (1 == idx)
		% mean sea level
		bc(1,1,idx).var = 'z';
		bc(1,1,idx).rhs = 0;
		% Dirichlet condition
		bc(1,1,idx).p   = 1;
		% wave entering from left
		bc(1,2,idx).var = 'z';
		bc(1,2,idx).rhs = z10;
		bc(1,2,idx).p   = [1,0];
		bc(1,2,idx).q   = [1,0];
	else
		bc(1,1,idx).var = '';
		bc(1,1,idx).rhs = [];
		bc(1,2,idx).var = '';
		bc(1,2,idx).rhs = [];
	end

	% upstream boundary condition
	if (1 == idx)
		bc(2,1,idx).var = '';
		bc(2,1,idx).rhs = [];
		bc(2,2,idx).var = '';
		bc(2,2,idx).rhs = [];
	else
		% river discharge
		bc(2,1,idx).var = 'Q';
		bc(2,1,idx).rhs = Q0;
		bc(2,1,idx).p   = 1;
		% wave entering from right / reflected wave
		bc(2,2,idx).var = 'z';
		bc(2,2,idx).rhs =   0;
		bc(2,2,idx).p   = [1,0];
		bc(2,2,idx).q   = [0,1];
	end
end % for idx

	% TODO, solve and compare to single channel

	hydrosolver    = BVPS_Characteristic();
	hydrosolver.xi = Xi;
	hydrosolver.nx = nx;
	opt.xs         = opt.xs;
	opt.dischargeisvariable = true;

	rtn = River_Tide_BVP( ...
				   'zb',      zb ...
				 , 'cd',      Cd ...
				 , 'width',   w0 ...
				 , 'omega',       omega ...
				 , 'opt',         opt ...
				 , 'hydrosolver', hydrosolver ...
				);
	rtn.bc = bc;

function [cid,eid,p] = jfun()
	cid = [1,2];
	eid = [2,1];
	p   = [1./w0];
end

	rtn.junction_condition = {@jfun}

	rtn.init();

	rtn.solve();

	figure(1)
	clf();
	for idx=1:rtn.nc
		x = rtn.x(idx);

		%subplot(length(rtn.rt),4,4*(idx-1)+1);
		subplot(2,4,1)
		z0 = rtn.z(0,idx);
		plot(x,z0);
		hold on
		
		subplot(2,4,5)
		w0 = rtn.width(idx);
		plot(x,w0);
		hold on
		
		%subplot(length(rtn.rt),4,4*(idx-1)+2);
		subplot(2,4,2)
		Q0 = rtn.Q(0,idx);
		plot(x,Q0);
		hold on
	
		[Qlr,zlr] = rtn.decompose(idx);
		k = 1;
		subplot(1,4,3)
		%subplot(length(rtn.rt),4,4*(idx-1)+3);
		z = rtn.z(k,idx);
		plot(x,abs(zlr(:,1)));
		%plot(x,abs(z));
		hold on
		set(gca,'colororderindex',1);
		plot(x,abs(zlr(:,2)),'--');
		set(gca,'colororderindex',1);
		plot(x,abs(z),':');
		
		subplot(1,4,4)
		%subplot(length(rtn.rt),4,4*(idx-1)+4);
		Q = rtn.Q(k,idx);
		plot(x,real(Q));
		hold on
		set(gca,'colororderindex',idx);
		plot(x,imag(Q),'--');
	end % for idx
end % test_river_tide_network_03

