% Fri 15 Jul 17:35:03 CEST 2016
% Karl Kastner, Berlin
%
%% Tide table
classdef Tidetable < handle
	properties
	% primary TPXO values
		time
		level
		ux
		uy
		umag
		udir
		dt

	% analysis values
		dmax
		dmin
		drange
		hdx
		ldx
		neap_dx
		range24
		spring_dx
		t24
		th
		tl
		tneap
		tspring
		vh
		vl
		x
		y
		placename
	end
	methods
		function obj = Tidetable()
		end
	end
end % class Tidetable

