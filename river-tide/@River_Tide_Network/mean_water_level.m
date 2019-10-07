% Fri 22 Feb 14:22:25 CET 2019
%% predict the mean water level
function [x,z0,Q0] = mean_water_level(obj,id)
		% point indices of channel end points
		channel = obj.channel;

		% point id of the selected channel
		pid = channel(id,:);

		% water level at terminal points and junctions
		z0    = obj.z0;

		% discharge through the channel
		Q0    = obj.Q0(id);

		% water level at channel end points
		z0  = z0(pid);
		
		L = obj.L(id);

		% grid for evalutation
		nx = max(2,round(L/obj.dx));
		x  = linspace(0,L,nx)';

		z0 = z0(1) + x/L*(z0(2)-z0(1));
		Q0 = Q0*ones(nx,1);
end % River_Tide_Network/mean_water_level

