% Sun  8 Oct 13:08:39 CEST 2017
% Karl Kastner, Berlin
%
% forcing by along-channel width-variation
function f = odefun_width(obj,f,k,w0,dw0_dx)
	% sign changed after test case 2
	f(:,2) = f(:,2) + -1./w0.^2.*dw0_dx;
end

