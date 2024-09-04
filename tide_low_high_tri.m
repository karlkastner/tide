% Mon  4 Nov 19:15:55 +08 2019
function [rr,y] = tide_low_high_tri(cc)
	c0     =  cc(:,1);
	c      =  cc(:,2:5);
	% ce = [0.*ce(:,1),1i*ce(:,2),-1i*ce(:,3),2i*ce(:,4),-2i*ce(:,5)]
	%cc = [0.*cc(:,1),cc(:,2:3),2*cc(:,4:5)]
	% to exp
	%ce = fourier_tri2exp(cc) 
	% derive
	%ce = [0.*ce(:,1),1i*ce(:,2),-1i*ce(:,3),2i*ce(:,4),-2i*ce(:,5)]
	% reorder
	p      =  [ +2i*c(:,4) - 2*c(:,3), ...
                    +1i*c(:,2) -   c(:,1), ...
                    zeros(length(c0),1), ...
                    +1i*c(:,2) +   c(:,1), ...
                    +2i*c(:,4) + 2*c(:,3)];
	% TODO avoid loops
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
			% TODO use roots4, there is a bug in roots4 when real part of z2 = 0
			z1_ = roots(p(idx,:));
		end
		z1(idx,:)  =  z1_;
	end
	r      = -1i*log(z1);
	fdx    =  real(r) < 0;
	r(fdx) =  r(fdx) + 2*pi;
	rr     =  real(r);
	ri     =  imag(r);

	y = zeros(n,4);
	% for each of the four zeros
	for idx=1:4
		%A = [cos(2*pi*rr), + sin(2*pi*rr), + cos(4*pi*rr), + sin(4*pi*rr)];
		%y      =  c0 + A*c.';
		y(:,idx) = c0 + cos(rr(:,idx)).*c(:,1) + sin(rr(:,idx)).*c(:,2) ...
			      + cos(2*rr(:,idx)).*c(:,3) + sin(2*rr(:,idx)).*c(:,4);
	end
end


