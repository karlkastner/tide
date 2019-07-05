% Sun  8 Oct 12:35:23 CEST 2017
%
%% tidal component of surface elevation determined from tidal discharge
%%
%% by continuity
%%
%% %   dz/dt + dq/dx = 0
%% => i o z = - dq/dx
%% =>     z = -1/(io) dq/dx
%% =>     z = 1i/o dq/dx
%%
%% TODO allow Q as input
function z = q_to_z(x,q,omega)
	z = 1i./omega*derivative1(x,q);
end

