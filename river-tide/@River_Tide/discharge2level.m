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
%% TODO rename into Qt_to_zt
%% Mon  7 Oct 19:04:14 PST 2019 : added correction for change of width
function z = discharge2level(obj,Q,x)
	if (nargin()<3)
		x     = obj.x;
	end
	w     = obj.width(x);
	omega = obj.omega;
	%q = bsxfun(@rdivide,Q,w);
	z = zeros(size(Q));
	D1_dx = obj.D1_dx(x);
	for idx=1:size(Q,2)
		z(:,idx) = 1i./(idx*omega*w) .* (D1_dx*Q(:,idx));
		%derivative1(x,Q(:,idx));
	    	  %*(         derivative1(x,q) ...
	          % + q./w.*derivative1(x,w) );
	end
end % River_Tide/q1_to_z1

