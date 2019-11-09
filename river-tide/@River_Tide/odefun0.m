% Sun 11 Mar 13:17:36 CET 2018
%% coefficients of the backwater equation for the river tide
%% TODO merge with backwater
function c = odefun0(obj,x,z0,z1,zb,w,dw_dx,Q0,Qhr,Qt,cd,fc)
	if (nargin() < 2)
		c = zeros(0,0,1);
		return;
	end

		g = obj.g;
		h    = z0-zb;
		% TODO no magic numbers
		hmin = 1;
		%h = max(h,hmin); %sqrt(eps));
		% sign of friction term is +, because minus sign is in fc(2)
		%h = (up(h) + 2*h + down(h))/4;
		A  = w.*h;
		Qt2 = sum(abs(Qt).^2,2);
		c        =  [zeros(size(x)), ...
			     g*A, ...
			     zeros(size(x)), .... 
			     ( -g*h.^2.*(h.^2 + 1/2*abs(z1).^2).*dw_dx ...
			       +cd.*w./(pi*A.^2).*(  ...
					      fc(:,1).*Qhr.^2 ... 
					    + fc(:,2).*Q0.*Qhr ...
					    + fc(:,3).*(abs(Q0).*Q0 + 1/2*Qt2) )) ...
			     ];

		% TODO
		% This quick fix recovering zero depth does not work well,
		% spurious oscillations emerge
if (1)
		fdx      = h<0;
		fdx      = find(fdx);
		c(fdx,:) = 0;
		%c(fdx,1) = 1;
		c(fdx,3) = 1;
		r = min(length(x),fdx+1);
		l = min(1,fdx-1);
		dzb_dx = (zb(r)-zb(l))./(x(r)-x(l));
		if (isscalar(w))
			w = repmat(w,length(x),1);
		end
		if (isscalar(cd))
			cd = repmat(cd,length(x),1);
		end
		h = (sqrt(cd(fdx)./g).*abs(Q0)./(w(fdx).*sqrt(abs(dzb_dx)))).^(2/3);
		c(fdx,4) = -(zb(fdx)+h);
end

end % River_Tide/odefun0

