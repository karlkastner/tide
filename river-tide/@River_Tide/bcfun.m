% Sun  8 Oct 13:08:39 CEST 2017
%% Robin (mixed) boundary conditions for the river tide,
%% supplied for each frequency component,
%% wrapper that copies values are from the member struct "bc"
%%
%%       q*(p*Q_1^- + (1-p)*dQ_1^-/dx 
%  + (1-q)*(p*Q_1^+ + (1-p)*dQ_1^+/dx) = rhs
%% input :
%%	x       : coordinate (left or right end)
%%	id,ccdx : frequency component index
%%                (1 = 0 omega (mean), 2 : 1 omega, 3 : 2 omega, ...
%% columns of bc : frequency
%% rows of bc, left, right boundary
%% output :
%% 	p : [2x1] linear combination of Dirichlet and Neumann boundary condition
%%          p(1) -> weight Dirichlet boundary condition
%%	    p(2) -> weight Neumann boundary condition
%%      q linear combination of left and right travelling (incoming and outgoing) wave
%%	    q(1) weight left going wave
%%          q(2) weight right going wave
%% 	rhs = 0 -> homogeneous boundary condition
%%
%%
%% function [rhs, p, q, obj] = bcfun(obj,x,y,ccdx)
function [rhs, p, q, obj] = bcfun(obj,x,y,ccdx)
	Xi = obj.Xi;
	switch (obj.opt.hmode)
	case {'matrix'}
		[rhs, p, q] = bc(x,y,ccdx);
	case {'iterate'}
		[rhs, p, q] = bc(x,y,ccdx+1);
	otherwise
		error('bcfun');
	end % switch

	function [rhs, p, q] = bc(x,~,id)
		switch (x)
		% TODO end should better be selected by index 1,2, rather than xl ,xr
		case {Xi(1)}
			rhs = obj.bc(1,id).rhs;
			p   = obj.bc(1,id).p;
			q   = obj.bc(1,id).q;
		case {Xi(2)}
			rhs = obj.bc(2,id).rhs;
			p   = obj.bc(2,id).p;
			q   = obj.bc(2,id).q;
		otherwise
			error('bcfun');
		end % switch
	end % bc
end % River_tide/bcfun

