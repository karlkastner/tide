% Fri 14 Aug 15:48:04 +08 2020
% coefficients of the backwater equation
% for coupling z0 and Q0 for determining the latter
function c = odefunQ0(obj,x,h0,z1,zb,w,dw_dx,Q0,Qhr,Qt,Cd,fc)
	g   = obj.g;
	A   = w.*h0;
	Qt2 = sum(abs(Qt).^2,2);
	fcw = +Cd.*w./(pi*A.^2).*fc;
	c = [  g*A, ...                                     % z'
	      zeros(size(x)), ...			    % z
	    -(  fcw(:,3).*abs(Q0) + fcw(:,2).*Qhr), ...       % Q0
	    -(  fcw(:,1).*Qhr.^2  + fcw(:,3).*0.5.*Qt2) ...  % inhomogeneous part
            ];
end % River_Tide/odefunQ0

