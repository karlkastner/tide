% Tue 19 Jul 15:46:05 CEST 2016
% Karl Kastner, Berlin
%% time of high water and low water
%% c.f. savenije 2012
function [eps_l, eps_hw] = savenije_timing_hw_lw()
	% eq 2.120
	eps_hw = atan( (1-delta*gamma) * (lambda/gamma + F*(phi*pi - 2*delta/gamma) ) );
	% eq 2.122
	eps_lw = atan( (1-delta*gamma) * (lambda/gamma + 2*F*phi*pi) );
end

