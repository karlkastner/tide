% Mon  4 Nov 19:15:55 +08 2019
function [r,y] = tide_low_high_exp(c_exp)
	if (size(c_exp,2) == 3)
		c_tri     = [c_exp(:,1),real(c_exp(:,2)),-imag(c_exp(:,2)),real(c_exp(:,3)),-imag(c_exp(:,3))];
	else
		c_tri     = [c_exp(:,1),real(c_exp(:,2)),-imag(c_exp(:,2))];
		% TODO, quick hack
		c_tri(:,4:5) = 0;
	end
	[r,y] = tide_low_high_tri(c_tri);
%	c_d = [0,conj(c_exp(2)),conj(2*c_exp(3))];
%	[r0,ry] = tide_slack_exp(c_d)
end

