% Sun  8 Oct 12:35:23 CEST 2017
% Karl Kastner, Berlin
%
%% determines tidal water surface amplitude (non-zero freqyency components of surface elevation)
%% from tidal discharge (non-zero freqyency components of the discharge)
%%
%% by continuity :
%%
%% dz/dt + dq/dx = 0
%% => i o z = - dq/dx
%% =>     z = -1/(io) dq/dx
%% =>     z = 1i/o dq/dx
%%
function z = discharge2level(obj,x,Qt,w)
	z = zeros(size(Qt));
	% TODO construct from basis functions, rather than derivative1
	dQt_dx = derivative1(x,Qt);
	for k=1:size(Qt,2)
		omega_k = obj.omega(k);
		z(:,k) = 1i./(omega_k*w) .* dQt_dx(:,k);
	end
end % River_Tide_BVP/q1_to_z1

