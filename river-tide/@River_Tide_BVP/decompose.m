% Fri 15 Dec 16:58:52 CET 2017
% Karl Kastner, Berlin
%
%% decompose the tide into a right and left travelling wave,
%% i.e. into incoming and reflected wave
%
% function [Q1lr, z1lr, obj] = decompose(obj)
%
%% TODO subtract forcing term
function [Q1lr, z1lr, obj] = decompose(obj,cdx)
	if (nargin()<2)
		cdx = 1;
	end
	z0    = obj.z(0,cdx);
	z1    = obj.z(1,cdx);
	Q0    = obj.Q(0,cdx);
	Q1    = obj.Q(1,cdx);
	x     = obj.x(cdx);
	w0    = obj.width(cdx,x);
	omega = obj.omega;

	[Q1lr, z1lr] = decompose@River_Tide(obj,x,w0,z0,z1,Q0,Q1);
end % decompose

