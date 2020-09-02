% Wed  2 Sep 10:37:55 +08 2020
function [out] = test_river_tide_network_02(rt_map,pflag)
	if (nargin()<2)
		pflag = 1;
	end
	out.id  = 2;
	out.name = 'single channel, split in two (two channels connected at ends)'

	% tidal amplitude
	z10  = 1;

	% river discharge
	Q0   = -10;

	% width of channel
	w0 = 1;

	% drag/friction coefficient
	Cd        = 2.5e-3;

	% depth of channel
	h0        = 10;

	% slope of channel bed
	S0        = -normal_flow_slope(Q0,h0,w0,drag2chezy(Cd));

	% bed level of channel
	zb     = @(x) -h0 + S0*x;

	% base frequency
	T         = Constant.SECONDS_PER_DAY;
	omega     = 2*pi/T;
	% domain size
	L = h0./S0;
	Xi        = L*[0,0.5;
                       0.5,1];

	% number of grid points
	nx               = [50,50];

	meta             = river_tide_test_metadata();
	opt              = meta.opt;
	opt.maxiter       = 100;
	opt.sopt.maxiter = 100;

	bc        = struct();
	for idx=1:2

	% downstream boundary condition
	if (1 == idx)
		% mean sea level
		bc(1,1,idx).var = 'z';
		bc(1,1,idx).rhs = 0;
	else
		bc(1,1,idx).var = '';
		bc(1,1,idx).rhs = [];
	end
	% Dirichlet condition
	bc(1,1,idx).p   = 1;

	% upstream boundary condition
	if (1 == idx)
		bc(2,1,idx).var = '';
		bc(2,1,idx).rhs = [];
	else
		% river discharge
		bc(2,1,idx).var = 'Q';
		bc(2,1,idx).rhs = Q0;
	end
	bc(2,1,idx).p   = 1;

	% wave entering from left
	if (1 == idx)
		bc(1,2,idx).var = 'z';
		bc(1,2,idx).rhs = z10;
	else
		bc(1,2,idx).var = '';
		bc(1,2,idx).rhs = [];
	end
	bc(1,2,idx).p   = [1,0];
	bc(1,2,idx).q   = [1,0];

	% wave entering from right / reflected wave
	if (1 == idx)
		bc(2,2,idx).var = '';
		bc(2,2,idx).rhs = [];
	else
		bc(2,2,idx).var = 'z';
		bc(2,2,idx).rhs =  0;
	end
	bc(2,2,idx).p   = [1,0];
	bc(2,2,idx).q   = [0,1];

	end % for idx

	% solve individually
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
	p   = [1./w0,1./w0];
end
	rtn.junction_condition = {@jfun}

	rtn.init();
	rtn.solve();

	figure(1)
	clf();
	for idx=1:rtn.nc
		x = rtn.x(idx);

		%subplot(length(rtn.rt),4,4*(idx-1)+1);
		subplot(1,4,1)
		z0 = rtn.z(0,idx);
		plot(x,z0);
		hold on
		
		%subplot(length(rtn.rt),4,4*(idx-1)+2);
		subplot(1,4,2)
		Q0 = rtn.Q(0,idx);
		plot(x,Q0);
		hold on
		k = 1;
		subplot(1,4,3)
		%subplot(length(rtn.rt),4,4*(idx-1)+3);
		z = rtn.z(k,idx);
		plot(x,real(z));
		hold on
		set(gca,'colororderindex',idx);
		plot(x,imag(z),'--');
		
		subplot(1,4,4)
		%subplot(length(rtn.rt),4,4*(idx-1)+4);
		Q = rtn.Q(k,idx);
		plot(x,real(Q));
		hold on
		set(gca,'colororderindex',idx);
		plot(x,imag(Q),'--')
	end % for idx

end % river_tide_test_06

