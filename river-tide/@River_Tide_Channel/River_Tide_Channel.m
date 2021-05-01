% Wed 11 Oct 10:56:52 CEST 2017
classdef River_Tide_Channel < River_Tide
	properties
		% index of the channel
		id
		hydrosolver
		rt_bvp;

		% water level, one column per frequency
		z
		zc
 		% discharge
		Q
		Qc
		% net sediment transport
		Qs

		% hydrodynamic boundary conditions
		% TODO, this should also go into "channel"
		bc = [struct('p',       1, 'rhs',  0, 'q', []   ,'var','z'), ... % 0 mean (wl or discharge) left
		      struct('p',       1, 'rhs',  1, 'q', [1 1],'var','Q'), ... % 1 main species left
		      struct('p', [0 1 0], 'rhs',  0, 'q', [1 1],'var','Q'), ... % 2 even overtide left
		      struct('p', [0 1 0], 'rhs',  0, 'q', [1 1],'var','Q'), ... % 3 triple overtide left
		      struct('p', [0 1 0], 'rhs',  0, 'q', [1 1],'var','Q'); ... % 4 quadruple overtide left
		      struct('p', [1 0 0], 'rhs', [], 'q', []   ,'var','Q'), ... % 0 mean (wl or discharge) right
		      struct('p', [1 0 0], 'rhs',  0, 'q', [1 1],'var','Q'), ... % 1 main species right
		      struct('p', [1 0 0], 'rhs',  0, 'q', [1 1],'var','Q'), ... % 2 even overtide right
		      struct('p', [1 0 0], 'rhs',  0, 'q', [1 1],'var','Q'), ... % 3 overtide right
		      struct('p', [1 0 0], 'rhs',  0, 'q', [1 1],'var','Q'), ... % 4 overtide right
		];

		% functions for channel properties
		fun = struct(  'z0', [] ...
			     , 'zb', [] ...
			     , 'width', [] ...
			     , 'cd', [] ...
			    );

		% temporary storage
		tmp = struct( 'D1_dx', [] ...
			      , 'D2_dx2', [] ...
			      ,'zb',[] ...
			      ,'width',[] ...
			      ,'c', struct( 'D1_dx',[] ...
			                   ,'D2_dx2',[] ...
					   ,'zb',[] ...
					   ,'width',[]) ...
			    );
	end % properties

	methods

	function obj = River_Tide_Channel(varargin)
                for idx=1:2:length(varargin)
			switch(varargin{idx})
			case {'fun.zb','zb','bed-level'}
				obj.set_zb(varargin{idx+1});
			case {'fun.width','width','w'}
				obj.set_width(varargin{idx+1});
			case {'fun.cd','cd'}
				obj.set_cd(varargin{idx+1});
			otherwise
                            obj = setfield_deep(obj,varargin{idx},varargin{idx+1});
			end % switch
		end % for idx
	end % constructor

	function nx = nx(obj)
		%nx = length(obj.hydrosolver.out(obj.id).x);
		nx = obj.hydrosolver.nx(obj.id);
		%out(obj.id).x);
	end

	function neq = neq(obj)
		neq = obj.hydrosolver.neq;
	end

	% TODO, the hydrosolver should neither store x, not out
	function x = x(obj)
%		if (nargin()<2)
%			cid = 1;
%		end
		x = obj.hydrosolver.out(obj.id).x;
	end
	
	function Xi = Xi(obj,eid)
		Xi = obj.hydrosolver.xi(obj.id,:);
		if (nargin>1)
			Xi = Xi(eid);
		end
	end

	function set_cd(obj,fun)
		if (isnumeric(fun))
			% constant value
			obj.fun.cd = @(x,~) fun*ones(size(x));
		elseif (isstr(fun))
			% lookup table
			tab = readtable(fun);
			x   = tab.x;
			cd  = tab.cd;
			fun = @(xi) interp1(x,cd,xi,'linear');
		elseif(isa(fun,'function_handle'))
		    switch (nargin(fun))
			case {1}
				obj.fun.cd = @(x,~) feval(fun,x);
    			case {2}
				obj.fun.cd = @(x,h) feval(fun,x,h);
%    			case {3}
%				obj.fun.cd = @(cid,x,h) feval(fun,x,h);
    		        otherwise
			    error('Drag coefficient function must take 1 (x) or 3 (x and h) arguments.');
		    end % switch nargin
	        else
 		    error('Drag coefficient must be a scalar, a function or a csv-filename');
		end % else of isa function
	end % set_cd

	function set_zb(obj,fun)
	    if (isnumeric(fun))
		obj.fun.zb = @(x) fun*ones(size(x));
	    elseif(isstr(fun))
		% lookup table
		tab = readtable(fun);
		x   = tab.x;
		zb  = tab.cd;
		fun = @(xi) interp1(x,zb,xi,'linear');
	    elseif(isa(fun,'function_handle'))
			obj.fun.zb = @(x) fun(x);
	    else
 		    error('Bed-Level must be a scalar, a function or a csv-filename');
	    end % else of isa function
	end % set_zb

	function set_width(obj,fun)
	    if (isnumeric(fun))
		obj.fun.width = @(x) fun*ones(size(x));
	    elseif(isstr(fun))
		% lookup table
		tab = readtable(fun);
		x   = tab.x;
		zb  = tab.width;
		fun = @(xi) interp1(x,width,xi,'linear');
	    elseif(isa(fun,'function_handle'))
			obj.fun.width = @(x) fun(x);
	   else
 		error('Bed-Level must be a scalar, a function or a csv-filename');
	   end % of else
	end % set_width

%	function out = out(obj)
%		out = obj.hydrosolver.out;
%	end

	function y = h0(obj, varargin)
%		if (nargin() < 3)
%			x = obj.x(cid);
%		end
		y = obj.waterlevel(0,varargin{:}) - obj.zb(varargin{:});
		%if (nargin()>1)
		%	y = interp1(obj.x,y,x,'linear');
		%end
	end

	function y = waterlevel(obj,fid,x)
%		if (nargin() < 3)
%			cid = 1;
%		end
		y = obj.z(:,fid+1);
		if (nargin()>2)
			y = interp1(obj.x,y,x,'linear');
		end
	end

	function y = A(obj,fid,varargin)
%		if (nargin() < 3)
%			cid = 1;
%		end
		w  = obj.width(varargin{:});
		z  = obj.waterlevel(fid,varargin{:}); %channel(cid).z;
		%(cid,varargin{:});
		zb = obj.zb(varargin{:});

		% TODO only when fid == 0
		z(:,1) = z(:,1) - zb;
		y  = w.*z;

		%y = obj.out(cid).z(:,fid+1);
		if (nargin()>2)
			y = interp1(obj.x,y,x,'linear');
		end
	end

	function y = discharge(obj,fid,x)
%		if (nargin()<3)
%			cid = 1;
%		end
		if (nargin()<2 || isempty(fid))
			y = obj.Q;
		else
			% TODO, problem
			y = obj.Q(:,fid+1);
		end	
		if (nargin()>3)
			y = interp1(obj.x,y,x,'linear');
		end
	end

	% discharge per unit width
	function y = q(obj,fid,x)
		if (nargin()<4)
			x = obj.x;
		end
		w = obj.width(x);
		y = obj.discharge(fid)./w;
		if (nargin()>2)
			y = interp1(obj.x,y,x,'linear');
		end
	end

	% TODO not correct with more than one frequency component
	function y = zrange(obj,x)
		y = 2*abs(obj.waterlevel(1));
		if (nargin()>1)
			y = interp1(obj.x,y,x,'linear');
		end
	end
	% TODO not correct with more than one frequency component
	function y = Qrange(obj,x)
		y = 2*abs(obj.discharge(1));
		if (nargin()>1)
			y = interp1(obj.x,y,x,'linear');
		end
	end
	% TODO not correct with more than one frequency component
	function y = zmid(obj,x)
		y = obj.waterlevel(0);
		if (nargin()>1)
			y = interp1(obj.x,y,x,'linear');
		end
	end
	% TODO not correct with more than one frequency component
	function y = Qmid(obj,x)
		y = obj.discharge(0);
		if (nargin()>1)
			y = interp1(obj.x,y,x,'linear');
		end
	end

	% velocity
	function y = u(obj,fid,varargin)
%		if (nargin() < 3)
%			cid = 1;
%		end
		y  = obj.q(fid,varargin{:})./obj.h0(varargin{:});
		%if (nargin()>2)
		%	y = interp1(obj.x,y,x,'linear');
		%end
	end

	% water level time series for one day
	function y = zt(obj,t)
%		if (nargin()<3)
%			cid = 1;
%		end
		y = obj.waterlevel(0);
		for k=1:obj.neq-1
			y = ( y + obj.waterlevel(k)*exp(2i*pi*k*obj.rt_bvp.omega*t) ...
			        + conj(obj.waterlevel(k)).*exp(-2i*pi*k*obj.omega*t) );
		end
		y = real(y);
	end

	% velocity time series for one day
	function y = ut(obj,t)
%		if (nargin()<3)
%			cid = 1;
%		end
		y = obj.u(0);
		omega = obj.rt_bvp.omega;
		for k=1:obj.neq-1
			y = ( y + obj.u(k)*exp(2i*pi*k*omega*t) ...
			        + conj(obj.u(k)).*exp(-2i*pi*k*omega*t) );
		end
		y = real(y);
	end % ut

	function y = velocity(obj,varargin)
		y = obj.u(varargin{:});
	end

	% tidal energy of first frequency component
	% TODO this should go into the River_Tide class
	% simply integrate over hole tidal cycle?
	function E = energy(obj,varargin)
		h0 = obj.depth(varargin{:});
		w0 = obj.width(varargin{:});
		z1 = abs(obj.waterlevel(1,varargin{:}));
		u1 = abs(obj.u(1,varargin{:}));
		E = energy@River_Tide(h0,w0,z1,u1);
	end % energy

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

	% cd may depend on the depth and thus cannot be precomputed
	function cd = cd(obj,x,h0)
%		if (nargin()<2)
%			cid = 1;
%		end
		if (nargin()<2)
			x = obj.x;
		end
		if (nargin()<3)
			h0 = obj.h0(x);
		end
		cd = obj.fun.cd(x,h0);
	end

	function zb = zb(obj,x)
%		if (nargin()<2)
%			cid=1;
%		end
		if (nargin()<2)
			nx = obj.nx;
		else
			nx = length(x);
		end
		switch (nx)
			case {obj.nx} %hydrosolver.nx(cid)}
				if (isempty(obj.tmp.zb))
					obj.tmp.zb = obj.fun.zb(x);
				end
				zb = obj.tmp.zb;
			case {obj.nx-1} %obj.hydrosolver.nx(cid)-1}
				if (isempty(obj.tmp.c.zb))
					obj.tmp.c.zb = obj.fun.zb(x);
				end
				zb = obj.tmp.c.zb;
			otherwise
				zb = obj.fun.zb(x);
		end % switch
	end % zb

	function width = width(obj,x)
%		if (nargin()<2)
%			cid = 1;
%		end
		if (nargin()<2)
			x = obj.x;
		end
		switch (length(x))
			case {obj.nx}
				if (isempty(obj.tmp.width))
					obj.tmp.width = obj.fun.width(x);
				end
				width = obj.tmp.width;
			case {obj.nx-1}
				if (isempty(obj.tmp.c.width))
					obj.tmp.c.width = obj.fun.width(x);
				end
				width = obj.tmp.c.width;
			otherwise
				width = obj.fun.width(x);
		end % switch
	end % width

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

	function D1_dx = D1_dx(obj,x)
		if (nargin()<2)
			x = obj.x;
		end
		switch (length(x))
		case {obj.nx}
			if (isempty(obj.tmp.D1_dx))
				obj.tmp.D1_dx = derivative_matrix_1_1d(x,[],2);
			end
			D1_dx = obj.tmp.D1_dx;
		case {obj.nx-1}
			if (isempty(obj.tmp.c.D1_dx))
				obj.tmp.c.D1_dx = derivative_matrix_1_1d(x,[],2);
			end
			D1_dx = obj.tmp.c.D1_dx;
		otherwise	
			D1_dx = 0;
		end
	end % D1_dx

	function D2_dx2 = D2_dx2(obj,x)
		if (nargin()<1)
			x = obj.x;
		end
		switch (length(x))
		case {obj.nx}%hydrosolver.nx(cid)}
			if (isempty(obj.tmp.D2_dx2))
				obj.tmp.D2_dx2 = derivative_matrix_2_1d(x); %,[],2);
			end
			D2_dx2 = obj.tmp.D2_dx2;
		case {obj.nx-1}
			if (isempty(obj.tmp.c.D2_dx2))
				obj.tmp.c.D2_dx2 = derivative_matrix_2_1d(x);%,[],2);
			end
			D2_dx2 = obj.tmp.c.D2_dx2;
		otherwise	
			D2_dx2 = 0;
		end
	end % D2_dx2

	function clear(obj)
		% to avoid assinging by dissimilar-structure
		obj.tmp = struct( 'D1_dx',[] ...
					   ,'D2_dx2', [] ...
			      ,'zb',[] ...
			      ,'width',[] ...
			      ,'c', struct( 'D1_dx', [] ...
					   ,'D2_dx2', [] ...
					   ,'zb',    [] ...
					   ,'width', [] ...
					  ) ...
			    );
		f_C = {'z','zc','Q','Qc','Qs'}
		for idx = 1:length(f_C)
			obj.(f_C{idx}) = [];
		end
	end % function clear
	end % methods
end % classded River_Tide_Channel

