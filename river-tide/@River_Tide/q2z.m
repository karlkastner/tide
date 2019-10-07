% Sun  8 Oct 12:35:23 CEST 2017
%
%% tidal component of surface elevation determined from tidal discharge
%%
%% by continuity :
%%
%% dz/dt + dq/dx = 0
%% => i o z = - dq/dx
%% =>     z = -1/(io) dq/dx
%% =>     z = 1i/o dq/dx
%%
%% TODO allow Q as input
%% TODO rename into Q1_to_z1
%% Mon  7 Oct 19:04:14 PST 2019 : added correction for change of width
function z = q_to_z(x,q,omega)
	w = obj.width(x);
	z = 1i./omega ...
	    *(         derivative1(x,q) ...
	       + q./w.*derivative1(x,w) );
end % River_Tide/q1_to_z1

