% Sun  8 Oct 13:08:39 CEST 2017
% Karl Kastner, Berlin
%
%% forcing by along-channel width-variation
% note : sign inverted to eq 3.4 in Kastner 2019
function f = odefun_width(obj,f,k,w0,dw0_dx)
	% Q'
	f(:,2) = f(:,2) + -1./(w0.*w0).*dw0_dx;
end

