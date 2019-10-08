% Fri  5 May 16:16:56 CEST 2017
% Karl Kastner, Berlin
%% process tidal data to extrac the tidal envelope
% TODO, this is partially a duplicate to Tide_Table
%
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
classdef Tidal_Envelope < handle
	properties
		% centre time for each day
		tc
		% high water and high water slack
		zh
		tzh
		tsh

		% low water and low water slack
		zl
		tzl
		tsl

		% highest velocity
		uh
		tuh

		% lowest velocity
		ul
		tul

		t0    = 0;
		Ti    = 25/24;
		order = 2;
	end % properties
	methods
		function obj = Tidal_Envelope(varargin)
			for idx=1:2:length(varargin)-1
				field = varargin{idx};
				val   = varargin{idx+1};
				obj.(field) = val;	
			end % for
		end % constructor
		
	function [zm obj] = zmid(obj)
		zm = 0.5*(obj.zh + obj.zl);
	end
	function [zr obj] = zrange(obj)
		zr = obj.zh - obj.zl;
	end
	function [um obj] = umid(obj)
		um = 0.5*(obj.uh + obj.ul);
	end
	function [ur obj] = urange(obj)
		ur = obj.uh - obj.ul;
	end
	function [dt obj] = dt_hws(obj)
		dt = obj.tsl - obj.tzh;
	end
	function [dt obj] = dt_lws(obj)
		dt = obj.tsl - obj.tzl;
	end
	function dphi_hws = dphi_hws(obj)
		% phase difference
		dphi_hws = 2*pi*(obj.dt_hws);
		dphi_hws = wrapToPi(dphi_hws);
	end
	function dphi_lws = dphi_lws(obj)
		dphi_lws = 2*pi*(dt_lws);
		dphi_lws = wrapToPi(dphi_lws);
	end
	
	end % methods
end % class

