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
function z = discharge2level(obj,x,Qt,w)
%	if (nargin()<3)
%		x     = obj.x(cdx);
%	end
%	w     = obj.width(cdx,x);
	omega = obj.omega;
	%q = bsxfun(@rdivide,Q,w);
	z = zeros(size(Qt));
%	D1_dx = obj.D1_dx(cdx,x);
	dQt_dx = derivative1(x,Qt);
	for idx=1:size(Qt,2)
		z(:,idx) = 1i./(idx*omega*w) .* dQt_dx(:,idx);
		%derivative1(x,Q(:,idx));
	    	  %*(         derivative1(x,q) ...
	          % + q./w.*derivative1(x,w) );
	end
end % River_Tide_BVP/q1_to_z1

