% Wed 11 Oct 10:56:52 CEST 2017
% Karl Kastner, Berlin
%
%% physical functions for computation of river tides in a single 1D channel
%% combined with BVP-solver in child-classes to determine the hydrodynamics
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
		pi = pi;

		% angular frequency of the main tidal species
		omega

		% options
		% TODO dischargeisvariables is options for bvp-solver
		opt = struct( ...
			       'friction_method', 'dronkers' ...
			     , 'friction_order', 2 ...
			     , 'oflag', [true,false(1,3)] ...
			     , 'dischargeisvariable', false ...
			    );
%		neq = 0;

		% TODO make part of opt
		flag = struct('aa',0,'gh',0,'oh',0);

		% symbolic computation
		issym           = false;

	end % properties
	methods
	    function obj = River_Tide(varargin)
                for idx=1:2:length(varargin)
			switch(varargin{idx})
			case {'opt'}
			    % this keeps default options that are not set
		            obj.opt = copy_fields(varargin{idx+1},obj.opt);
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
	    end % River_Tide (constructor)

	end % methods
end % class River_Tide

