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
	x     = obj.x;
	z0    = obj.waterlevel(0);
%	z1    = obj.waterlevel(1);
	zb    = obj.zb();
%	Q0    = obj.discharge(0);
%	Q1    = obj.discharge(1);
	w0    = obj.width();
	Cd    = obj.cd(z0-zb);

	%[Q1lr, z1lr] = obj.rt.decompose(x,w0,z0,z1,zb,Q0,Q1,Cd);
	[Q1lr, z1lr] = obj.rt.decompose(x,obj.Q,obj.Q,zb,w0,Cd);
end % decompose

