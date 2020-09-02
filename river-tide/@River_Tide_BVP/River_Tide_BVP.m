% Wed 11 Oct 10:56:52 CEST 2017
%
%% hydrodynamics and morphodynamics of 1D tidal channel networks
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
% TODO integrate seasonal variation into boundary conditions
% TODO implement multigrid-like solver for morphodynamics
% TODO move non-linear forcing terms to lhs, let odefun return the coefficients
classdef River_Tide_BVP < River_Tide
	properties
		% water level and discharge, one column per frequency
		out = struct('z',[],'zc',[],'Q',[],'Qc',[]);

		% hydrodynamic solver
		hydrosolver;

		% hydrodynamic boundary conditions
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

		% boundary condition for sediment transport (kg/s)
		bc_Qs;

		% functions for channel properties
		fun = struct(  'z0', [] ...
			     , 'zb', [] ...
			     , 'width', [] ...
			     , 'cd', [] ...
			    );
		
		% minimum water depth (m)
		hmin = 0.1;

		% sediment properties for transport
		sediment = struct( 'd_mm', 0.2, ... % m      grain diameter
				   'p',    0.6, ... % 1      packing density
				   'rho', 2650 ...  % kg/m^3 material density
				 );

		% time stepper for determining morphodynamics
		morsolver;

		% object performing transport coupling at junctions
		bifurcation;

		% hydrodynamic coupling conditions at channel junctions
		junction_condition = {};

		% sediment transport coupling conditions
		junction_Qs

		% history over time
		evolution = struct('t', [], 'zb', []);
		
		% temporary storage
		tmp = struct( 'D1_dx',[] ...
			      ,'zb',[] ...
			      ,'width',[] ...
			      ,'c', struct( 'D1_dx',[] ...
					   ,'zb',[] ...
					   ,'width',[]) ...
			    );
		cflag;
	end % properties
	methods
	function obj = River_Tide_BVP(varargin)
		obj.opt.imethod      = 'spline';
		obj.opt.iorder       = 1;
		obj.opt.stokes_order = 2;
		obj.bifurcation      = Bifurcation();
		obj.bifurcation.division_rule = @obj.bifurcation.sediment_division_geometric;
		obj.hydrosolver      = BVPS_Characteristic();

                for idx=1:2:length(varargin)
			switch(varargin{idx})
			case {'hydrosolver'}
				obj.hydrosolver = varargin{idx+1};
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
			case {'fun.zb','zb','bed-level'}
				obj.set_zb(varargin{idx+1});
			case {'fun.width','width','w0'}
				obj.set_width(varargin{idx+1});
			case {'fun.cd','cd'}
				obj.set_cd(varargin{idx+1});
			otherwise
                            obj = setfield_deep(obj,varargin{idx},varargin{idx+1});
			end % switch
                end %for idx
		obj.hydrosolver.inifun = @obj.initial_value;
	end % River_Tide (constructor)

	function set_cd(obj,fun)
	if(isa(fun,'function_handle'))
		switch (nargin(fun))
			case {1}
				obj.fun.cd = @(~,x,h) feval(fun,x);
    			case {2}
				obj.fun.cd = @(cid,x,h) feval(fun,x);
    			case {3}
				obj.fun.cd = @(cid,x,h) feval(fun,x,h);
    		otherwise
			error('Chezy coefficient function must take 3 (cid,x) or 3 (cid, x and h) arguments.');
		end
	else
		% constant value
		%obj.fun.cd = @(cdx,x,~) fun(cdx)*ones(size(x));
		if (isscalar(fun))
			obj.fun.cd = @(~,x,~) fun*ones(size(x));
		else
			obj.fun.cd = @(cdx,x,~) fun(cdx)*ones(size(x));
		end
	end % else of isa function
	end
	function set_zb(obj,fun)
	if(isa(fun,'function_handle'))
		switch (nargin(fun))
		case {1}
			obj.fun.zb = @(~,x) fun(x);
		case {2}
			obj.fun.zb = @(cdx,x) fun(cdx,x);
		otherwise
			error('bed-level must be a vector with one entry per channel or a function accepting 1 (x) or 2 (cid,x) argumets');
		end
	else
		obj.fun.zb = @(cdx,x) fun(cdx)*ones(size(x));
	end
	end


	function set_width(obj,fun)
	if(isa(fun,'function_handle'))
		switch (nargin(fun))
		case {1}
			obj.fun.width = @(~,x) fun(x);
		case {2}
			obj.fun.width = @(cdx,x) fun(cdx,x);
		otherwise
			error('width must be a vector with one entry per channel or a function accepting 1 (x) or 2 (cid,x) argumets');
		end
	else
		if (isscalar(fun))
			obj.fun.width = @(~,x) fun*ones(size(x));
		else
			obj.fun.width = @(cdx,x) fun(cdx)*ones(size(x));
		end
	end
	end
	% quasi/pseudo members
	function nc = nc(obj)
		nc = obj.hydrosolver.nc;
	end

	function nx = nx(obj)
		nx = obj.hydrosolver.nx;
	end

	function neq = neq(obj)
		neq = obj.hydrosolver.neq;
	end

	function x = x(obj,cid)
		if (nargin()<2)
			cid = 1;
		end
		x = obj.hydrosolver.out(cid).x;
	end
	
	function Xi = Xi(obj,cid,eid)
		Xi = obj.hydrosolver.xi(cid,:);
		%Xi = obj.opt.Xi(cid, :);
		if (nargin>2)
			Xi = Xi(eid);
		end
	end

%	function out = out(obj)
%		out = obj.hydrosolver.out;
%	end

	function y = h0(obj, cid, x)
		if (nargin() < 3)
			x = obj.x(cid);
		end
		y = obj.z(0,cid,x) - obj.zb(cid, x);
		%if (nargin()>1)
		%	y = interp1(obj.x,y,x,'linear');
		%end
	end

	function y = z(obj,fid,cid,x)
		if (nargin() < 3)
			cid = 1;
		end
		y = obj.out(cid).z(:,fid+1);
		if (nargin()>3)
			y = interp1(obj.x(cid),y,x,'linear');
		end
	end
	
	function y = Q(obj,fid,cid,x)
		if (nargin()<3)
			cid = 1;
		end
		y = obj.out(cid).Q(:,fid+1);
		if (nargin()>3)
			y = interp1(obj.x,y,x,'linear');
		end
	end

	function y = q(obj,cid,fid,x)
		if (nargin()<4)
			x = obj.x;
		end
		w = obj.fun.width(cid,x);
		y = obj.Q(cid,fid)./w;
		if (nargin()>2)
			y = interp1(obj.x,y,x,'linear');
		end
	end


	% TODO consider q1 not work any more with two frequency components
	function y = zrange(obj,cid,x)
		y = 2*abs(obj.z(1,cid));
		if (nargin()>1)
			y = interp1(obj.x,y,x,'linear');
		end
	end
	function y = Qrange(obj,fid,cid,x)
		y = 2*abs(obj.Q(1,cid));
		if (nargin()>1)
			y = interp1(obj.x,y,x,'linear');
		end
	end
	function y = zmid(obj,cid,x)
		y = obj.z(0,cid);
		if (nargin()>1)
			y = interp1(obj.x,y,x,'linear');
		end
	end
	function y = Qmid(obj,cid,x)
		y = obj.Q(0,cid);
		if (nargin()>1)
			y = interp1(obj.x,y,x,'linear');
		end
	end

	function y = u(obj,fid,cid,varargin)
		if (nargin() < 3)
			cid = 1;
		end
		y  = obj.q(cid,fid,varargin{:})./obj.h0(cid,varargin{:});
		%if (nargin()>2)
		%	y = interp1(obj.x,y,x,'linear');
		%end
	end

	function y = velocity(obj,varargin)
		y = obj.u(varargin{:});
%		Q = obj.Q_;
%		if (nargin()>1)
%			Q = interp1(obj.x,Q,x,'linear');
%		else
%			x = obj.x;
%		end
%		w = obj.fun.width(x);
%		h0 = obj.h0(x);
%		y = bsxfun(@times,Q,1./(w.*h0));
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
	function cd = cd(obj,cid,x,h0)
		if (nargin()<3)
			x = obj.x(cid);
		end
		if (nargin()<4)
			h0 = obj.h0(cid,x);
		end
		cd = obj.fun.cd(cid,x,h0);
	end

	function zb = zb(obj,cid,x)
		if (nargin()<2)
			cid=1;
		end
		if (nargin()<3)
			nx = obj.nx;
		else
			nx = length(x);
		end
		switch (nx)
			case {obj.hydrosolver.nx(cid)}
				if (isempty(obj.tmp(cid).zb))
					obj.tmp(cid).zb = obj.fun.zb(cid,x);
				end
				zb = obj.tmp(cid).zb;
			case {obj.hydrosolver.nx(cid)-1}
				if (isempty(obj.tmp(cid).c.zb))
					obj.tmp(cid).c.zb = obj.fun.zb(cid,x);
				end
				zb = obj.tmp(cid).c.zb;
			otherwise
				zb = obj.fun.zb(cid,x);
		end
	end

	function width = width(obj,cid,x)
		if (nargin()<2)
			cid = 1;
		end
		if (nargin()<3)
			x = obj.x(cid);
		end
		switch (length(x))
			case {obj.hydrosolver.nx(cid)}
				if (isempty(obj.tmp(cid).width))
					obj.tmp(cid).width = obj.fun.width(cid,x);
				end
				width = obj.tmp(cid).width;
			case {obj.hydrosolver.nx(cid)-1}
				if (isempty(obj.tmp(cid).c.width))
					obj.tmp(cid).c.width = obj.fun.width(cid,x);
				end
				width = obj.tmp(cid).c.width;
			otherwise
				width = obj.fun.width(cid,x);
		end
	end % width

	function D1_dx = D1_dx(obj,cid,x)
		if (nargin()<2)
			x = obj.x;
		end
		switch (length(x))
		case {obj.hydrosolver.nx(cid)}
			if (isempty(obj.tmp(cid).D1_dx))
				obj.tmp(cid).D1_dx = derivative_matrix_1_1d(x,[],2);
			end
			D1_dx = obj.tmp(cid).D1_dx;
		case {obj.hydrosolver.nx(cid)-1}
			if (isempty(obj.tmp(cid).c.D1_dx))
				obj.tmp(cid).c.D1_dx = derivative_matrix_1_1d(x,[],2);
			end
			D1_dx = obj.tmp(cid).c.D1_dx;
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
		% to avoid assinging by dissim-structure
		obj.tmp = struct( 'D1_dx',[] ...
			      ,'zb',[] ...
			      ,'width',[] ...
			      ,'c', struct( 'D1_dx',[] ...
					   ,'zb',[] ...
					   ,'width',[]) ...
			    );
		for idx=2:obj.nc
		obj.tmp(idx) = struct( 'D1_dx',[] ...
			      ,'zb',[] ...
			      ,'width',[] ...
			      ,'c', struct( 'D1_dx',[] ...
					   ,'zb',[] ...
					   ,'width',[]) ...
			    );
		end % for idx
		obj.out = struct();
	end % clear

	end % methods
end % class River_Tide_BVP

