% Tue  5 Nov 09:42:23 +08 2019
function [rr,y,r] = tide_slack_exp(c)
	% real(c0  + c1 e1 + c2 e2)
	% c0  + 1/2(c1 e1 + c1' e-1) + 1/2 (c2 e2 + c' e-2) = 0
	% 2 c0 e2 + c1 e3 + c1' e1 + c2 e4 + c2' e0 = 0
	p = [c(:,3),c(:,2),2*c(:,1),conj(c(:,2)),conj(c(:,3))];
	n = size(c,1);
	z = zeros(n,4);
	for idx=1:n
		z(idx,:) = rvec(roots(p(idx,:)));
	end
	r = -1i*log(z)/(2*pi);
	r = real(r);
	fdx = r<0;
	r(fdx) = r(fdx)+1;
	fdx = r>1;
	r(fdx) = r(fdx)-1;
	rr = r;
	y = c(:,1)*ones(1,4);
	for idx=1:2
		y = y + bsxfun(@times,c(:,idx+1),exp(idx*2i*pi*rr));
	end
end

