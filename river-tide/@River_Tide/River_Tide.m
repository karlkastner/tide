% Wed 11 Oct 10:56:52 CEST 2017
%% river tide in a single 1D channel
%% TODO split in two classes:
%% one that stores data (RT_Solve), one that provides equations (RT_Analytic)
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
classdef River_Tide < handle
	properties
		% acceleration by gravity
		g = Constant.gravity;

		% bed level
		zb
		% width

		%wfun
		w
		% drag coefficient
		%cdfun

		% angular frequency of tide
		omega

		% domain start and end point
		Xi

		% discretisation of Xi
		x
		% TODO store Q and z as matrices
		% water level, one column per frequency
		z_
		%zmid
		%zrange
		% discharges, one column per frequency
		%Q0fun
		Q_
		%qmid
		%qrange

		% Backwater1D object to solve for tidally averaged water level
		backwater

		% numerical options
		opt = struct( 'nx',     1024 ...
			     ,'xs',     10 ... % stretch factor of mesh
			     ,'model_str',  'wave' ...
			     ,'solver', @bvp2fdm ...
			     ... %,'solver', @bvp2fdm ...
			     ,'friction_method', 'dronkers' ...
			     ,'friction_order', 2 ...
			     ,'hmode', 'iterate' ... % 'no-tide' ...
			     , 'imethod', 'spline' ...
			     , 'o2', false ...
			    );

		% downstream boundary condition
		z0_downstream = 0; % [0 sqrt(eps)];
		% river discharge
		Q0_
		bc = [struct('p', 1, 'rhs',  0, 'q', []   ,'var','z'), ... % mean (wl or discharge) left
		      struct('p', 1, 'rhs',  1, 'q', [1 1],'var','Q'), ... % main species left
		      struct('p', [0 1 0], 'rhs',  0, 'q', [1 1],'var','Q'); ... % even overtide ;
		      struct('p', [1 0 0], 'rhs', [], 'q', []   ,'var','Q'), ... % mean (wl or discharge) right
		      struct('p', [1 0 0], 'rhs',  0, 'q', [1 1],'var','Q'), ... % main species right
		      struct('p', [1 0 0], 'rhs',  0, 'q', [1 1],'var','Q')  ... % even overtide right
			 ];

		flag = struct('aa',0,'gh',0,'oh',0);

		% convergence flag
		cflag
		
		tmp = struct( 'D1',[]);
		fun = struct( 'z0',[] ...
			     ,'zb',[] ...
			     ,'width',[] ...
			     ,'cd',[]);

		% finite-volume object,
		% if the tide is computed by solving the full SWE
		fv
		% quantities in time, if full SWE
		T
		H
		U
		
		initial;

		issym           = false;
	end % properties
	methods (Static)
		%z1 = q2z(x,q1,omega1);
		[f, F3]  = odefun1(Q0, Qhr, Q1, Q2, h0, dh0_dx, dz0_dx, w, dw_dx, cd, g, c, omega1, flag);
		[f ]     = odefun2(Q0, Qhr, Q1, Q2, h0, dh0_dx, dz0_dx, w, dw_dx, cd, g, c, omega1, flag);
%		[k0, k] = wave_number_analytic(Q0,W,H,cd,omega,az1,Qt);
	end % static
	methods
	function obj = River_Tide(varargin)
                for idx=1:2:length(varargin)
			switch(varargin{idx})
			case {'opt'}
			    % this keeps default options that are not set
		            obj.opt = copy_fields(varargin{idx+1},obj.opt);
			case {'bc'}
				bc = varargin{idx+1};
				for jdx=1:length(bc)
				    obj.bc(jdx) = copy_fields(bc(jdx),obj.bc(jdx));
				end
			case {'issym'}
				obj.issym = varargin{idx+1};
				if (obj.issym)
					syms g positive
					obj.g = g;
				end
			otherwise
                            obj = setfield_deep(obj,varargin{idx},varargin{idx+1});
			end
                end %for idx
		% object has to be created here,
		% as otherwise matlab does not allow for parallel objects
		if (isempty(obj.backwater))
			obj.backwater = Backwater1D('rt',obj);
		end
	end % River_Tide (constructor)

	% quasi members
	function y = h0(obj, x)
		y = obj.z(0)-obj.zb;
		if (nargin()>1)
			y = interp1(obj.x,y,x);
		end
	end

	function y = z(obj,id,x)
		y = obj.z_(:,id+1);
		if (nargin()>2)
			y = interp1(obj.x,y,x);
		end
	end
	
	function y = Q(obj,id,x)
		y = obj.Q_(:,id+1);
		if (nargin()>2)
			y = interp1(obj.x,y,x);
		end
	end

	function y = q(obj,id,x)
		y = obj.Q(id)./obj.w;
		if (nargin()>2)
			y = interp1(obj.x,y,x);
		end
	end
	function y = u(obj,id,x)
		y  = obj.q(id)./obj.h0;
		if (nargin()>2)
			y = interp1(obj.x,y,x);
		end
	end

%	function y = z0(obj,x)
%		y = obj.z(:,1);
%		if (nargin()>1)
%			y = interp1(obj.x,y,x);
%		end
%	end
%	function y = z1(obj,x)
%		y = obj.z(:,2);
%		if (nargin()>1)
%			y = interp1(obj.x,y,x);
%		end
%	end
%	function y = z2(obj,x)
%		y = obj.z(:,3);
%		if (nargin()>1)
%			y = interp1(obj.x,y,x);
%		end
%	end
%
%	function y = Q0(obj,x)
%		y = obj.Q(:,1);
%		if (nargin()>1)
%			y = interp1(obj.x,y,x);
%		end
%	end
%
%	function y = Q1(obj,x)
%		y = obj.Q(:,2);
%		if (nargin()>1)
%			y = interp1(obj.x,y,x);
%		end
%	end
%
%	function y = Q2(obj,x)
%		y = obj.Q2(:,3);
%		if (nargin()>1)
%			y = interp1(obj.x,y,x);
%		end
%	end
%
%	function y = q0(obj, x)
%		y = obj.Q0./obj.w;
%		if (nargin()>1)
%			y = interp1(obj.x,y,x);
%		end
%	end
%
%	function y = q1(obj, x)
%		y = obj.Q1./obj.w;
%		if (nargin()>1)
%			y = interp1(obj.x,y,x);
%		end
%	end
%
%	function y = q2(obj, x)
%		y = obj.Q2./obj.w;
%		if (nargin()>1)
%			y = interp1(obj.x,y,x);
%		end
%	end
%	
%	function y = u0(obj, x)
%		y  = obj.q0./obj.h0;
%		if (nargin()>1)
%			y = interp1(obj.x,y,x);
%		end
%	end
%
%	function [y obj] = u1(obj, x)
%		y   = obj.q1./obj.h0;
%		if (nargin()>1)
%			y = interp1(obj.x,y,x);
%		end
%	end
%
%	function [y obj] = u2(obj, x)
%		y   = obj.q2./obj.h0;
%		if (nargin()>1)
%			y = interp1(obj.x,y,x);
%		end
%	end

	% TODO consider q1 not work any more with two frequency components
	function y = zrange(obj,x)
		y = 2*abs(obj.z1);
		if (nargin()>1)
			y = interp1(obj.x,y,x);
		end
	end
	function y = Qrange(obj,x)
		y = 2*abs(obj.Q1);
		if (nargin()>1)
			y = interp1(obj.x,y,x);
		end
	end
	function y = zmid(obj,x)
		y = obj.z0;
		if (nargin()>1)
			y = interp1(obj.x,y,x);
		end
	end
	function y = Qmid(obj,x)
		y = obj.Q0;
		if (nargin()>1)
			y = interp1(obj.x,y,x);
		end
	end
	function y = velocity(obj,x)
		Q = obj.Q_;
		if (nargin()>1)
			Q = interp1(obj.x,Q,x);
		else
			x = obj.x;
		end
		w = obj.fun.width(x);
		h0 = obj.tmp.h0(x);
		y = bsxfun(@times,Q,1./(w.*h0));
	end

	function [dy_dx, obj] = dy_dx(obj,field,x)
		if (iscell(field))
			if (length(field)>1)
				y = obj.(field{1})(field(2:end));
			else
				y = obj.(field{1});
			end
		else
			y     = obj.(field);
		end
		%dy_dx = cdiff(y)./cdiff(obj.x);
		dy_dx = derivative1(obj.x,y);
		if (nargin()>2)
			dy_dx = interp1(obj.x,dy_dx,x);
		end
	end

	function [absy, obj] = amplitude(obj, field, id, x)
		if (iscell(field))
			if (length(field)>1)
				y = obj.(field{1})(field(2:end));
			else
				y = obj.(field{1});
			end
		else
			y     = obj.(field)(id);
		end
		absy  = abs(y);
		if (nargin()>3)
			absy = interp1(obj.x,absy,x);
		end
	end
	
	function [zs,obj] = animate(obj,T)
		zs = repmat(obj.z0,1,length(T));
		for idx=1:length(T)
			zs(:,idx) = real(obj.z1.*exp(1i*obj.omega*T(idx)));
		end
	end

	function [ay, obj] = admittance(obj, varargin)
		absy = obj.amplitude(varargin{:});
		ay = absy./absy(1);
	end

	function [k, obj] = wave_number(obj, field, id)
		% z = hat z exp(i(ot-kx)) => dz/dx = -i k z => k = 1i/z dz/dx
		% dQ/dx = -ikQ
		y     = obj.(field)(id);
		dy_dx = derivative1(obj.x,y);
		k     = 1i*(dy_dx./y);
	end

	function [py, obj] = phase(obj, field, id, x)
		if (iscell(field))
			if (length(field)>1)
				y = obj.(field{1})(field(2:end));
			else
				y = obj.(field{1});
			end
		else
			y     = obj.(field)(id);
		end
		py  = angle(y);
		if (nargin()>3)
			py = interp1(obj.x,py,x);
		end
	end
	function [r, obj] = damping_modulus(obj,x)
		az1 = obj.amplitude('z1');
		daz1_dx = obj.dy_dx({'amplitude','z1'});
		r   = daz1_dx./az1;
                %dz1_dx0(kdx,jdx) = (az1(2)-az1(1))/(x(2)-x(1))*1/arg.z1_downstream;
		if (nargin()>1)
			r = interp1(obj.x,r,x);
		end
	end

	function dz0_dx = dz0_dx(obj,x)
		x_      = obj.tmp.x;
		z0_     = obj.tmp.z0(x_);
		xc       = 0.5*(x_(2:end)+x_(1:end-1));
		dx       = x_(2:end)-x_(1:end-1);
		dz0      = z0_(2:end)-z0_(1:end-1);
		dz0_dx_c = dz0./dx;
		dz0_dx   = interp1(xc,dz0_dx_c,x,'pspline','extrap');
	end
	
	function dh0_dx = dh0_dx(obj,x)
		x_       = obj.tmp.x;
		h0_      = obj.tmp.h0(x_);

		xc       = 0.5*(x_(2:end)+x_(1:end-1));
		dx       = x_(2:end)-x_(1:end-1);
		dh0      = h0_(2:end)-h0_(1:end-1);
		dh0_dx_c = dh0./dx;
		dh0_dx  = interp1(xc,dh0_dx_c,x,'pspline','extrap');
	end

	end % methods
end % class River_Tide

