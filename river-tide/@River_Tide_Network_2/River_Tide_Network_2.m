% Mon 10 Aug 12:38:17 +08 2020
classdef River_Tide_Network_2 < handle
	properties
		rt
		nx
		m
	end
	methods
		function obj = River_Tide_Network_2(rt)
			if (nargin()>0)
				obj.rt = rt;
			end
		end

		% pseudo members
		function Q = Q(obj,id)
			Q = [];
			for idx=1:length(obj.rt)
				Q = [Q; obj.rt(idx).Q(id)];
			end	
		end
		
		function z = z(obj,id)
			z = [];
			for idx=1:length(obj.rt)
				z = [z; obj.rt(idx).z(id)];
			end	
		end	

		function init(obj)
			for idx=1:length(obj.rt)
				obj.rt(idx).init();
			end
		end
	end % methods
end % class River_Tide_Network

