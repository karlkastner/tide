% Sun  8 Oct 13:08:39 CEST 2017
% Karl Kastner, Berlin
%
%% Robin (mixed) boundary conditions for the river tide,
%% supplied for each frequency component,
%% wrapper that copies values are from the member struct "bc"
%%
%%       q*(p*Q_1^- + (1-p)*dQ_1^-/dx 
%  + (1-q)*(p*Q_1^+ + (1-p)*dQ_1^+/dx) = rhs
%% input :
%%	cid : channel index
%%	bif : 1,2 : index for letft/right end of channel
%%	fid : frequency component index
%%                (1 = 0 omega (mean), 2 : 1 omega, 3 : 2 omega, ... )
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
%% function [rhs, p, q, obj] = bcfun(obj,cid,bid,fid)
%
%function [rhs, p, q, set, obj] = bcfun(obj,cid,bid,fid,varargin)
function [rhs, p, q, type, obj] = bcfun(obj,cid,bid,fid,varargin)
	if (isa(obj.bc(bid,fid,cid).rhs,'function_handle'))
		rhs = feval(obj.bc(bid,fid,cid).rhs,varargin{:});
		% TODO, retransform here, if necessary
	else
		rhs = obj.bc(bid,fid,cid).rhs;
	end
	p    = obj.bc(bid,fid,cid).p;
	q    = obj.bc(bid,fid,cid).q;
	type = obj.bc(bid,fid,cid).var;
end % River_tide/bcfun

