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
% TODO opt.oflag at the moment cannot skip frequencies
%classdef River_Tide_Network < handle
classdef River_Tide_Network < handle
	properties
		% river tide object
		rt

		% channel data
		channel

		% hydrodynamic solver
		hydrosolver;


		% boundary condition for sediment transport (kg/s)
		% TODO, this should also go into "channel"
		bc_Qs;
		

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
		%evolution = struct('t', [], 'zb', []);

		% convergence flag		
		cflag;
	end % properties

	methods
	function obj = River_Tide_Network(varargin)
		% these default values can be overwritten with varargin
		% they cannot be set in "properties" as this leads
		% to concurrency problems between multiple instances


		obj.bifurcation      = Bifurcation();
		obj.bifurcation.division_rule = @obj.bifurcation.sediment_division_geometric;
		obj.hydrosolver      = BVPS_Characteristic();

		% note that inheriting does not work, due to call to rt functions by RT_Channel
		obj.rt = River_Tide();
		% cannot be set in properties, as opt is inherited
		obj.rt.opt.imethod = 'spline';
		obj.rt.opt.iorder  =  1;
		% expansion of stokes flow at bifurcations
		obj.rt.opt.stokes_order = 2;
		% change of grid spacing along channel
		obj.rt.opt.xs = 1;
		% minimum water depth (m)
		obj.rt.opt.hmin = 0.1;

                for idx=1:2:length(varargin)
			switch(varargin{idx})
			case {'hydrosolver'}
				obj.hydrosolver = varargin{idx+1};
			case {'opt'}
			    % this keeps default options that are not set
		            obj.rt.opt = copy_fields(varargin{idx+1},obj.rt.opt);
			case {'bc'}
				bc = varargin{idx+1};
				for jdx=1:length(bc)
				    obj.bc(jdx) = copy_fields(bc(jdx),obj.bc(jdx));
				end
			case {'omega'}
				obj.rt.omega = varargin{idx+1};
%			case {'issym'}
%				obj.issym = varargin{idx+1};
%				if (obj.issym)
%					syms g positive
%					obj.g = g;
%				end
			otherwise
                            obj = setfield_deep(obj,varargin{idx},varargin{idx+1});
			end % switch
                end %for idx


	end % River_Tide (constructor)

	function transform_bc(obj)
		for cdx=1:obj.nc
			obj.channel(cdx).transform_bc();
		end
	end

	function [rhs, p, q, type] = bcfun(obj,cid,varargin)
		[rhs, p, q, type] = obj.channel(cid).bcfun(varargin{:});
	end

	function y0 = initial_value(obj,cid,x)
		y0 = obj.channel(cid).initial_value(x);
	end

	function [f, obj] = odefun(obj,cdx,varargin)
		f = obj.channel(cdx).odefun(varargin{:});	
	end

	function [rmse,res] = check_continuity(obj)
		for cdx=1:length(obj.nc)
			[rmse(cdx),res{cdx}] = obj.channel(cdx).check_continuity();
		end
	end % River_Tide_Network / check_continuity

	% quasi/pseudo members
	function nc = nc(obj)
		%obj.hydrosolver.nc;
		nc = length(obj.channel);
	end

	%function nx = nx(obj)
	%	nx = obj.hydrosolver.nx;
	%end

	function neq = neq(obj)
		neq = obj.hydrosolver.neq;
	end

	function clear(obj)
		obj.cflag = [];
		for idx=1:obj.nc
			obj.channel(idx).clear(); % River_Tide_Channel.empty();
		end
		% to avoid assinging by dissimilar-structure
%		obj.tmp = struct( 'D1_dx',[] ...
%					   ,'D2_dx', [] ...
%			      ,'zb',[] ...
%			      ,'width',[] ...
%			      ,'c', struct( 'D1_dx', [] ...
%					   ,'D2_dx', [] ...
%					   ,'zb',    [] ...
%					   ,'width', [] ...
%					  ) ...
%			    );
%		for idx=2:obj.nc
%		obj.tmp(idx) = struct( 'D1_dx',[] ...
%					   ,'D2_dx', [] ...
%			      ,'zb',[] ...
%			      ,'width',[] ...
%			      ,'c', struct( 'D1_dx',[] ...
%					   ,'D2_dx',[] ...
%					   ,'zb',   [] ...
%					   ,'width',[] ...
%					  ) ...
%			    );
%		end % for idx
		% struct();
	end % clear

	end % methods
end % class River_Tide_Network

