% Thu 10 Oct 11:35:49 PST 2019
function [rmse,res] = check_continuity(obj)
	% TODO, derivative is directly available
	x = obj.x;
	omega = obj.omega;
	w = obj.fun.width(x);
	Q1 = obj.Q_(:,2);
	z1 = obj.z_(:,2);
	res = 1i*omega*w.*z1 + derivative1(x,Q1);
	rmse = rms(res);
end % check_continuity

