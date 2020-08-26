% Sun 11 Mar 13:17:36 CET 2018
%% coefficients of the backwater equation for the river tide
%% TODO merge with backwater
function c = odefunz0(obj,x,h0,z1,zb,w,dw_dx,Q0,Qhr,Qt,Cd,fc)
	if (nargin() < 2)
		c = zeros(0,0,1);
		return;
	end
	g   = obj.g;
	A   = w.*h0;
	Qt2 = sum(abs(Qt).^2,2);
	% note that there is no dw/dx term
	% sign of friction term is + because minus sign is in fc(2)


	fcw = +Cd.*w./(pi*A.^2);
	c       =  [   zeros(size(x)) ...
		     , g*A ...
		     , zeros(size(x)) .... 
			fcw.*(  ...			
				      fc(:,1).*Qhr.^2 ... 
				    + fc(:,2).*Q0.*Qhr ...
				    + fc(:,3).*(abs(Q0).*Q0 + 0.5*Qt2) ...
					  ) ...
		   ];
end % River_Tide/odefun0


