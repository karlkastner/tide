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
		% zb

		% angular frequency of the main tidal species
		omega

		% domain left and right end
		Xi

		% grid spanning dmoain Xi
		x

		% water level, one column per frequency
		z_
		zc_

		% discharges, one column per frequency
		Q_
		Qc_

		% river discharge (scalar)
		% TODO remove this
		Q0_

		% Backwater1D object to solve for tidally averaged water level
		% superfluous if matrix mode is chosen
		backwater

		% numerical options
		opt = struct( 'nx',     1024 ...
			     ,'xs',     10 ... % stretch factor of mesh
			     ,'model_str',  'wave' ...
			     ,'solver', @bvp2c ...
			     ... %,'solver', @bvp2fdm ...
			     ,'friction_method', 'dronkers' ...
			     ,'friction_order', 2 ...
			     ... %,'hmode', 'iterate' ... % 'no-tide' ...
			     , 'hmode', 'matrix' ... % 'no-tide' ...
			     , 'imethod', 'spline' ...
			     , 'oflag', [true,false(1,3)] ...
			     , 'ifun', [] ...
			     , 'dischargeisvariable', false ...
			    );

		% downstream boundary condition
		% TUDO remove, superfluous if matrix mode is chosen
		z0_downstream = 0;
		neq = 0;

		bc = [struct('p', 1, 'rhs',  0, 'q', []   ,'var','z'), ...       % 0 mean (wl or discharge) left
		      struct('p', 1, 'rhs',  1, 'q', [1 1],'var','Q'), ...       % 1 main species left
		      struct('p', [0 1 0], 'rhs',  0, 'q', [1 1],'var','Q'), ... % 2 even overtide left
		      struct('p', [0 1 0], 'rhs',  0, 'q', [1 1],'var','Q'), ... % 3 overtide left
		      struct('p', [0 1 0], 'rhs',  0, 'q', [1 1],'var','Q'); ... % 4 overtide left
		      struct('p', [1 0 0], 'rhs', [], 'q', []   ,'var','Q'), ... % 0 mean (wl or discharge) right
		      struct('p', [1 0 0], 'rhs',  0, 'q', [1 1],'var','Q'), ... % 1 main species right
		      struct('p', [1 0 0], 'rhs',  0, 'q', [1 1],'var','Q'), ... % 2 even overtide right
		      struct('p', [1 0 0], 'rhs',  0, 'q', [1 1],'var','Q'), ... % 3 overtide right
		      struct('p', [1 0 0], 'rhs',  0, 'q', [1 1],'var','Q'), ... % 4 overtide right
		];

		flag = struct('aa',0,'gh',0,'oh',0);

		% convergence flag
		out=struct();
		%cflag
		
		tmp = struct( 'D1_dx',[],'zb',[],'width',[], ...
			      'c', struct('D1_dx',[],'zb',[],'width',[]));

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

		hmin = 0.2;


		% sediment properties for transport
		sediment = struct( 'd_mm', 0.2, ... % m      grain diameter
				   'p',    0.6, ... % 1      packing density
				   'rho', 2650 ...  % kg/m^3 material density
				 );

		% time stepper for determining morphodynamics
		morsolver;

		% boundary condition for sediment transport (kg/s)
		bc_Qs;
	end % properties
	methods (Static)
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
			case {'fun.cd','cdfun','cd'}
				cdfun = varargin{idx+1};
				if(isa(cdfun,'function_handle'))
					switch (nargin(cdfun))
			    			case {1}
							obj.fun.cd = @(x,h) feval(cdfun,x);
			    			case {2}
							obj.fun.cd = @(x,h) feval(cdfun,x,h);
			    		otherwise
						error('Chezy coefficient function must take 1 (x) or 2 (x and h) arguments.');
					end
				else
					% constant value
					obj.fun.cd = @(x,h) cdfun;
				end % else of isa function
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
		if (nargin()<2)
			x= obj.x;
		end
		y = obj.z(0,x)-obj.zb(x);
		%if (nargin()>1)
		%	y = interp1(obj.x,y,x,'linear');
		%end
	end

	function y = z(obj,id,x)
		y = obj.z_(:,id+1);
		if (nargin()>2)
			y = interp1(obj.x,y,x,'linear');
		end
	end
	
	function y = Q(obj,id,x)
		y = obj.Q_(:,id+1);
		if (nargin()>2)
			y = interp1(obj.x,y,x,'linear');
		end
	end

	function y = q(obj,id,x)
		if (nargin()<3)
			x = obj.x;
		end
		w = obj.fun.width(x);
		y = obj.Q(id)./w;
		if (nargin()>2)
			y = interp1(obj.x,y,x,'linear');
		end
	end
	function y = u(obj,id,x)
		y  = obj.q(id)./obj.h0;
		if (nargin()>2)
			y = interp1(obj.x,y,x,'linear');
		end
	end


	% TODO consider q1 not work any more with two frequency components
	function y = zrange(obj,x)
		y = 2*abs(obj.z1);
		if (nargin()>1)
			y = interp1(obj.x,y,x,'linear');
		end
	end
	function y = Qrange(obj,x)
		y = 2*abs(obj.Q1);
		if (nargin()>1)
			y = interp1(obj.x,y,x,'linear');
		end
	end
	function y = zmid(obj,x)
		y = obj.z0;
		if (nargin()>1)
			y = interp1(obj.x,y,x,'linear');
		end
	end
	function y = Qmid(obj,x)
		y = obj.Q0;
		if (nargin()>1)
			y = interp1(obj.x,y,x,'linear');
		end
	end
	function y = velocity(obj,x)
		Q = obj.Q_;
		if (nargin()>1)
			Q = interp1(obj.x,Q,x,'linear');
		else
			x = obj.x;
		end
		w = obj.fun.width(x);
		h0 = obj.h0(x);
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
			dy_dx = interp1(obj.x,dy_dx,x,'linear');
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
			absy = interp1(obj.x,absy,x,'linear');
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
			py = interp1(obj.x,py,x,'linear');
		end
	end
	function [r, obj] = damping_modulus(obj,x)
		az1 = obj.amplitude('z1');
		daz1_dx = obj.dy_dx({'amplitude','z1'});
		r   = daz1_dx./az1;
                %dz1_dx0(kdx,jdx) = (az1(2)-az1(1))/(x(2)-x(1))*1/arg.z1_downstream;
		if (nargin()>1)
			r = interp1(obj.x,r,x,'linear');
		end
	end

	function dz0_dx = dz0_dx(obj,x)
		x_       = obj.x;
		z0_      = obj.z0(x_);
		xc       = 0.5*(x_(2:end)+x_(1:end-1));
		dx       = x_(2:end)-x_(1:end-1);
		dz0      = z0_(2:end)-z0_(1:end-1);
		dz0_dx_c = dz0./dx;
		dz0_dx   = interp1(xc,dz0_dx_c,x,'pspline','extrap');
	end
	
	function dh0_dx = dh0_dx(obj,x)
		x_       = obj.x;
		h0_      = obj.h0(x_);

		xc       = 0.5*(x_(2:end)+x_(1:end-1));
		dx       = x_(2:end)-x_(1:end-1);
		dh0      = h0_(2:end)-h0_(1:end-1);
		dh0_dx_c = dh0./dx;
		dh0_dx  = interp1(xc,dh0_dx_c,x,'pspline','extrap');
	end

	% cd may depend on the depth and thus cannot be precomputed
	function cd = cd(obj,x,h0)
		if (nargin()<2)
			x = obj.x;
		end
		if (nargin()<3)
			h0 = obj.h0(x);
		end
		cd = obj.fun.cd(x,h0);
	end

	function zb = zb(obj,x)
		switch (length(x))
			case {obj.opt.nx}
				if (isempty(obj.tmp.zb))
					obj.tmp.zb = obj.fun.zb(x);
				end
				zb = obj.tmp.zb;
			case {obj.opt.nx-1}
				if (isempty(obj.tmp.c.zb))
					obj.tmp.c.zb = obj.fun.zb(x);
				end
				zb = obj.tmp.c.zb;
			otherwise
				zb = obj.fun.zb(x);
		end
	end

	function width = width(obj,x)
		if (nargin()<2)
			x = obj.x;
		end
		switch (length(x))
			case {obj.opt.nx}
				if (isempty(obj.tmp.width))
					obj.tmp.width = obj.fun.width(x);
				end
				width = obj.tmp.width;
			case {obj.opt.nx-1}
				if (isempty(obj.tmp.c.width))
					obj.tmp.c.width = obj.fun.width(x);
				end
				width = obj.tmp.c.width;
			otherwise
				width = obj.fun.width(x);
		end
	end % width

	function D1_dx = D1_dx(obj,x)
		if (nargin()<2)
			x = obj.x;
		end
		switch (length(x))
		case {obj.opt.nx}
			if (isempty(obj.tmp.D1_dx))
				obj.tmp.D1_dx = derivative_matrix_1_1d(x,[],2);
			end
			D1_dx = obj.tmp.D1_dx;
		case {obj.opt.nx-1}
			if (isempty(obj.tmp.c.D1_dx))
				obj.tmp.c.D1_dx = derivative_matrix_1_1d(x,[],2);
			end
			D1_dx = obj.tmp.c.D1_dx;
		otherwise	
			D1_dx = 0;
		end
	end % D1_dx

	function [zs,obj] = animate(obj,T)
		zs = repmat(obj.z0,1,length(T));
		for idx=1:length(T)
			zs(:,idx) = real(obj.z1.*exp(1i*obj.omega*T(idx)));
		end
	end
	function clear(obj)
		obj.tmp = struct( 'D1_dx',[],'zb',[],'width',[], ...
			      'c', struct('D1_dx',[],'zb',[],'width',[]));
	end

	end % methods
end % class River_Tide

