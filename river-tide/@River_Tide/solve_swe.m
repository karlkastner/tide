% Fri 15 Dec 10:29:28 CET 2017
%% determine river tide by the fully non-stationary FVM and then extract the tide
%  this is experimental and not yet fuylly working
function solve_swe(obj)
			g      = 9.81;
			h0     = -zb_downstream;
			c      = sqrt(g*h0);
			% lambda = 2*pi*c/omega;
			Tp     = 2*pi/omega;
			Tx     = (Xi(2)-Xi(1))/c;

			Ti    = [0,2*Tx+Tp];
			%zbfun = @(x) -h0 + S0*x;
			zbfun_ = @(x) zbfun(x,zb_downstream,S0);
			wfun   = @(x) W0*ones(size(x));
			cdfun  = @(x) cd*ones(size(x));
			a      = z1_downstream;
			hin    = @(t) h0+a*sin(omega*t);
			bc{1}  = @(obj,t,q,dt) obj.bc_level(t, q, dt, hin, W0)
			% was - Q0
			bc{2}  = @(obj,t,q,dt) obj.bc_inflow(t, q, dt, @(t) -Q0, cd, W0);
			icfun  = 'backwater';
			% TODO, move to SWE

			% fully solve non-linear SWE
			[T Xi H U fv] = fv_swe(Ti, ...
					       Xi, ...
					       zbfun_, ...
					       wfun, ...
					       cdfun, ...
					       bc, ...
					       icfun, ...
					       Q0, ...
					       a, ...
					       opt);
			% post process
			% resample
			fdx = (T>=T(end)-Tp);
			T_  = T(fdx);
			n   = length(T);
			Tr  = linspace(T(1),T(end),n);
			dt  = (T(end)-T(1))/(n-1);
			H   = interp1(T,H',Tr)';
			U   = interp1(T,U',Tr)';

			% extract range over the last period		
			te  = Tidal_Envelope('Ti',86400,'order',0);
			
			zrange = zeros(size(H,1),1);
			zmid   = zeros(size(H,1),1);
			z0     = zeros(size(H,1),1);
			z1     = zeros(size(H,1),1);
			q0     = zeros(size(H,1),1);
			q1     = zeros(size(H,1),1);
			qrange = zeros(size(H,1),1);
			qmid   = zeros(size(H,1),1);
			for idx=1:size(H,1)
				te.init(Tr(1),dt,H(idx,:),H(idx,:).*U(idx,:)); %,U(idx,:));
				%te.init(Tr(1),dt,H(idx,:),U(idx,:));
				zrange_ = te.zrange;
				zmid_ = te.zmid;
				zrange(idx) = zrange_(end);
				zmid(idx)   = zmid_(end);
				qrange_     = te.urange();
				qrange(idx) = qrange_(end);
				qmid_       = te.umid();
				qmid(idx)   = qmid_(end);
			end
			

			% frequency components
			% TODO, maybe z and q are also better stored as vectors
			Ti_ = Ti/86400;
			stft = STFT('t0',Ti_(1),'tend',Ti_(end), ...
			             'Ti',1,'T',[1 1/2],'order',0);
			for idx=1:size(H,1)
			stft.transform(Tr/86400,H(idx,:));
			c  = stft.icoeff();
			% TODO end is nan
			z0(idx) = c(1,end-1);
			z1(idx) = c(2,end-1);

			stft.transform(Tr/86400,H(idx,:).*U(idx,:));
			c  = stft.icoeff();
			q0(idx) = c(1,end-1);
			q1(idx) = c(2,end-1);
			end

			% h to z
			zmid = zmid+fv.pde.zb;
			z0   = z0+fv.pde.zb;

			% check for convergence
			% difference between last two periods
			% STFT

			obj        = RT();
			obj.fv     = fv;
			obj.x      = fv.x;
			obj.z0     = z0;
			obj.z1     = z1;
			obj.zrange = zrange;
			obj.zmid   = zmid;
			obj.qrange = qrange;
			obj.qmid   = qmid;
			obj.Q0     = W0*q0;
			obj.q1     = q1;
			obj.zb     = fv.pde.zb;
			obj.w      = fv.pde.w;
			obj.cd      = fv.pde.cd;

			obj.T = T;
			obj.H = H;
			obj.U = U;
			% TODO
			% obj.zb    = fv.zb;
			% obj.z0    = z0;
			% obj.z1    = z1;
			%obj.q0
			%obj.q1
			%obj.
			%obj.cd
			%obj.w
end % solve_swe


