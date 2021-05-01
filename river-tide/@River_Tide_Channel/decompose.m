% Fri 15 Dec 16:58:52 CET 2017
% Karl Kastner, Berlin
%
%% decompose the tide into a right and left travelling wave,
%% i.e. into incoming and reflected wave
%
% function [Q1lr, z1lr, obj] = decompose(obj)
%
%% TODO subtract forcing term
function [Q1lr, z1lr, obj] = decompose(obj)
	if (nargin()<2)
		cdx = 1;
	end
	z0    = obj.waterlevel(0);
	z1    = obj.waterlevel(1);
	Q0    = obj.discharge(0);
	Q1    = obj.discharge(1);
	x     = obj.x;
	w0    = obj.width(x);
	omega = obj.rt_bvp.omega;

	[Q1lr, z1lr] = decompose@River_Tide(obj,x,w0,z0,z1,Q0,Q1);
end % decompose

