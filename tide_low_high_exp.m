% Mon  4 Nov 19:15:55 +08 2019
%
% min and max of 0.5*Re(c(2)*exp(2i*pi*t) + c(3)*exp(4i*pi*t))
%
function [r,y] = tide_low_high_exp(c_exp)
	if (isvector(c_exp))
		c_exp = rvec(c_exp);
	end
	%if (size(c_exp,2) == 3)
	%	c_tri     = [c_exp(:,1),2*real(c_exp(:,2)),-2*imag(c_exp(:,2)),2*real(c_exp(:,3)),-2*imag(c_exp(:,3))];
	%else
	%	c_tri     = [c_exp(:,1),2*real(c_exp(:,2)),-2*imag(c_exp(:,2))];
	%	% TODO, quick hack
	%	%c_tri(:,4:5) = 0;
	%end
	c_tri = fourier_exp2tri(c_exp);
	[r,y] = tide_low_high_tri(c_tri);
	y = real(y);
%	c_d = [0,conj(c_exp(2)),conj(2*c_exp(3))];
%	[r0,ry] = tide_slack_exp(c_d)
end

