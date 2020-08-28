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
%
%function [rhs, p, q, set, obj] = bcfun(obj,cid,bid,varargin)
function [rhs, p, q, set, obj] = bcfun(obj,cid,bid,varargin)
	if (isa(obj.bc(bid,cid).rhs,'function_handle'))
		rhs = feval(obj.bc(bid,cid).rhs,varargin{:});
		% TODO, retransform here, if necessary
	else
		rhs = obj.bc(bid,cid).rhs;
	end
	p   = obj.bc(bid,cid).p;
	q   = obj.bc(bid,cid).q;
	if (~obj.opt.dischargeisvariable)
		set = obj.bc(bid,cid).set;
	else
		set = obj.bc(bid,cid).var;
	end
end % River_tide/bcfun

