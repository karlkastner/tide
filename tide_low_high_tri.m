% Mon  4 Nov 19:15:55 +08 2019
function [rr,y] = tide_low_high_2(c)
	c0     =  c(:,1);
	c      =  c(:,2:5);
	p      =  [c(:,4)*2i - 2*c(:,3), + c(:,2)*1i - c(:,1), zeros(length(c0),1), + c(:,2)*1i + c(:,1), + c(:,4)*2i + 2*c(:,3)];
	% TODO avoid loops, use quartic function
	n    = size(c,1);
	z1   = NaN(n,4);
	for idx=1:n
		if (0 == p(idx,1))
			if (0 == p(idx,3))
				z1_ = NaN(1,4);
			else
				z1_ = roots(p(idx,3:end));
			end
		else
			z1_ = roots(p(idx,:));
		end
		z1(idx,:)  =  z1_;
	end
	r      = -log(z1)*1i;
	r      =  r/(2*pi);
	fdx    =  real(r) < 0;
	r(fdx) =  r(fdx)+1;
	rr     =  real(r);
	ri     =  imag(r);
%	rr = rr.'
%	rr = rr(:)

	y = zeros(n,4);
	% for each of the four zeros
	for idx=1:4
		%A = [cos(2*pi*rr), + sin(2*pi*rr), + cos(4*pi*rr), + sin(4*pi*rr)];
		%y      =  c0 + A*c.';
		y(:,idx) = c0 + cos(2*pi*rr(:,idx)).*c(:,1) + sin(2*pi*rr(:,idx)).*c(:,2) ...
			      + cos(4*pi*rr(:,idx)).*c(:,3) + sin(4*pi*rr(:,idx)).*c(:,4);
	end
end


