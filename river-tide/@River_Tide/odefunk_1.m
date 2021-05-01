% Sun  8 Oct 13:08:39 CEST 2017
% Karl Kastner, Berlin
% coefficients of the frequency components of the tide
% 
function [f]  = odefunk_1(obj, Q, QQ, Qhr, zs,zb, h0, dh0_dx, dz0_dx, w0, dw0_dx, Cd, cf, D1_dx, D2_dx)
	% start with index 2, index 1 : backwater is computed elsewhere
	k  = 2;
	for idx=1:length(obj.opt.oflag)
	    if (obj.opt.oflag(idx))
		f(:,:,k) = obj.odefunk_1_(idx, Q, QQ, Qhr, h0, dh0_dx, dz0_dx, ...
						    w0, dw0_dx, Cd, cf, D1_dx);
	    	k=k+1;
	    end % if
	end % for 
end % odefunk_1

