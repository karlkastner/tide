% Sun 20 Feb 17:30:13 CET 2022
% Note : the phase of the overtide is not correct without including the dQ/dt term!
% TODO	when we solve for zs, then we indirectly solve the full swe, stability? necessary to use wave at all?
% TODO check that junction condition is properly provided
classdef River_Tide_Wave_Equation < handle
	properties
		nx
		L
		w
		cd
		h
		g = 9.81;
		x
		D1 = {};
		D2 = {};
		dt
		nf = 4;
		T
		Q_ini
		bc
		flag = struct('aa', true);
		omega = 2*pi/86400;
		junction
		cid = [];
		% first and last index of channel
		first
		last
	end
	methods
		function obj = RTT(obj0)
		end

		function nc = nc(obj)
			nc = length(obj.nx);
		end

		function obj = init(obj)
			% TODO saninity checks of parameters
			obj.first = cumsum(obj.nx)-obj.nx+1;
			obj.last  = cumsum(obj.nx);
			obj.x = zeros(obj.last(end),1);
			for cdx=1:obj.nc
				obj.x(obj.first(cdx):obj.last(cdx)) = linspace(0,obj.L(cdx),obj.nx(cdx))';
				obj.D1{cdx} = derivative_matrix_1_1d(obj.nx(cdx),obj.L(cdx));
				obj.D2{cdx} = derivative_matrix_2_1d(obj.nx(cdx),obj.L(cdx));
			end
			obj.Q_ini = zeros(obj.last(end),1);
		end

		% first order coefficients for river tides
		%
		% 2Q/AdQ/dx - g h d^2Q/dx^2 + 2 cd w/A^2 abs(Q) dQ/dt = 0
		%
		% d^2Q/dt^2 ignored
		% advective acceleration neglected
		% rise of mean water level allong channels neglected
		function c = dQ_dt_c(rt,t,Q)
			c = zeros(obj.last(end),3);
			% for each channel
			for cdx=1:obj.nc
				c(obj.first(cdx):obj.last(cdx),:) = ...
                                     [ 2*rt.cd(cdx)./(rt.w(cdx).*rt.h(cdx).*rt.h(cdx)).*abs(Q(obj.first(cdx):obj.last(cdx))), ...
				       zeros(rt.nx(cdx),1), ...
                                      -rt.g*rt.h(cdx).*(rt.D2{cdx}*Q(obj.first(cdx):obj.last(cdx)))];
			end
			% note : the aa term cannot be included, as it contains a mixed derivative of form d/dx d/dt Q
			%if (rt.flag.aa)
			%	c(:,1) = c(:,1) + 2*Q./(rt.h.*rt.w).*(rt.D1*Q);
			%end
		end
	
		function dQ_dt = dQ_dt(rt,t,Q)
			c = rt.dQ_dt_c(t,Q);
			dQ_dt = -1./c(:,1).*(c(:,2).*Q + c(:,3));
		end

		% wave equation coefficients for river tides
		%
		% d^2Q/dt^2 - g A 1/w d^2Q/dx^2 + d/dt cd*w*Q|Q|/A^2
		% d^2Q/dt^2 - g h d^2Q/dx^2 + cd w/A^2 abs(Q) dQ/dt = 0
		%
		% advective acceleration optional
		% no width convergence
		% deep channel   : h constant time
		% no tidal flats : w constant in time
		function c = d2Q_dt2_c(rt,Q)
			c = [];
			for cdx=1:obj.nc
				c = [c; ones(rt.nx(cdx),1),
                                     2*rt.cd(cdx)./(rt.w(cdx).*rt.h(cdx).*rt.h(cdx)).*abs(Q(obj.first(cdx):obj.last(cdx))), ...
				     zeros(rt.nx(cdx),1), ...
				     -rt.g*rt.h(cdx).*(rt.D2{cdx}*Q(obj.first(cdx):obj.last(cdx)))];
			end % for cdx
		end % d2Q_dt2_c

		function dy_dt = d2Q_dt2(obj,t,y)
			Q       = y(1:obj.last(end));
			dQ_dt   = y(obj.last(end)+1:2*obj.last(end));
			d2Q_dt2 = [];
			for cdx=1:obj.nc
                        d2Q_dt2_ = -(-obj.g*obj.h(cdx).*(obj.D2{cdx}*Q(obj.first(cdx):obj.last(cdx))) ...
                                     + 2*obj.cd(cdx)./(obj.h(cdx).*obj.h(cdx).*obj.w(cdx)) ...
				       .*abs(Q(obj.first(cdx):obj.last(cdx))).*dQ_dt(obj.first(cdx):obj.last(cdx)));
			if (obj.flag.aa)
				d2Q_dt2_ = -2*(Q.*(obj.D1*dQ_dt) + dQ_dt.*(obj.D1*Q))./(obj.w.*obj.h);
				% TODO this is a quick fix as the bnd values are not properly defined
				%d2Q_dt2_(1) = 0;
				%d2Q_dt2_(end) = 0;
				d2Q_dt2 = d2Q_dt2 + d2Q_dt2_; 
			end % for cdx
				d2Q_dt2 = [d2Q_dt2;
                                           d2Q_dt2_];
			end % for cdx
			dy_dt = [dQ_dt;
                                 d2Q_dt2];
		end

		% transforms second order wave equation in set of two first order equations
		function [t,Q,dQ_dt,zs] = solve_d2Q_dt2(obj)
			y_ini = [obj.Q_ini; zeros(obj.last(end),1)];
			[t,y] = solve_dy_dt(@obj.d2Q_dt2,@obj.bcfun,obj.T,y_ini,obj.dt);
			Q     = y(1:obj.last(end),:);
			dQ_dt    = y(obj.last(end)+1:2*obj.last(end),:);
			if (nargout()>3)
				zs = obj.Q2zs(t,Q,dQ_dt);
			end
		end

		function zs = Q2zs(rt,t,Q,dQ_dt)
			if (rt.nc > 1)
				error('not yet implemented');
			end
			dQ_dt_ = diff(Q,[],2)./(t(2)-t(1));
			dQ_dt_ = inner2outer(dQ_dt_,2);
			% dQ/dt + aa + g h dzs_dx + cd w Q|Q|/A^2
			dzs_dx = mid(-(dQ_dt_ + rt.cd.*Q.*abs(Q)./(rt.w.*rt.h.*rt.h))./(rt.g*rt.h*rt.w));
			dx     = diff(rt.x);
			% TODO both sides
			cz0    = rt.bc(1).rhs;
			if (isa(cz0,'function_handle'))
				cz0 = cz0(t).';
			end
			zs0 = cz0(1,:);
			size(cz0)
			for idx=1:size(cz0,1)-1
				zs0 = (  zs0 ...
				       + cz0(idx+1,:).*exp(+1i*idx*rt.omega*rvec(t)) ...
				       + conj(cz0(idx+1,:)).*exp(-1i*idx*rt.omega*rvec(t)) ...
                                      );
			%sin(rt.omega*t);
			end
			zs     = cumsum([rvec(zs0);
                                         dzs_dx.*dx]);
		end

		function [t,Q] = solve_dQ_dt(obj)
			[t,Q] = solve_dy_dt(@obj.dQ_dt,@obj.bcfun,obj.T,obj.Q_ini,obj.dt);
			% y to z
			% dz_dt = dy_dx
			% dy/dt + g A dz_dx + cd w abs(y)|y|/A
			% z(:,tdx) =
		end

		function T = T1(rt)
			T = 2*pi/rt.omega;
		end
		function c = y2c(rt,t,y)
			A  = fourier_matrix_exp(1./(1:rt.nf)*rt.T1,t);
			c = (A \ y.').';
		end
		function [y,yt] = c2y(rt,t,c)
			A  = fourier_matrix_exp(1./(1:rt.nf)*rt.T1,t);
			y = real(A*c.').';
			yt = zeros([size(y),rt.nf+1]);
			if (nargout()>1)
			yt(:,:,1) = real(A(:,1)*c(:,1).').';
			for idx=1:rt.nf
				id = 2*idx;
				yt(:,:,idx+1) = real(A(:,id:id+1)*c(:,id:id+1).').';
			end
			end
		end

		function cz = cQ2cz(rt,cQ)
			if (rt.nc > 1)
				error('not yet implemented');
			end
			%A  = fourier_matrix_exp([1,1/2,1/3,1/4]*86400,t);
			cz = zeros(size(cQ,1)-1,size(cQ,2));
			for jdx=1:rt.nf
				% w dz/dt + dQ/dx = 0
				% dz/dt = -1/w dQ/dx
				% i j o (zj ep - conj(zj) em) = -1/w d(Qj ep + Qj em)/dx
				% (zj ep - conj(zj) em) = -1/(i o j w) d(Qj ep + conj(Qj) em)/dx
				% compare coeff :  
				% zj ep       = -1/(i o j w) d Qj ep
				% conj(zj) em = +1/(i o j w) d conj(Qj) em
				% zj = -1/(i o
				cz(:,2*jdx)   = -1./(+jdx*1i*rt.omega*rt.w)*diff(cQ(:,2*jdx))./diff(rt.x);
				cz(:,2*jdx+1) = -1./(-jdx*1i*rt.omega*rt.w)*diff(cQ(:,2*jdx+1))./diff(rt.x);
			end
		end

		% local channel index to global netwok index
		%function id = cid2nid(obj,cid,idx)
		%	sn = sum(obj.nx)-obj.nx;
		%	id = sn(cid)+idx;
		%end
	end % methods
end % RTT

