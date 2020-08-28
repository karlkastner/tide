% Wed 11 Oct 10:18:54 CEST 2017
% Karl Kastner, Berlin
%
%% extract values of individual variables from BVP-solver result vector
%
% TODO alternatively, keep all variables in a single vector and just once compute the indices
function [z0, Q0, Qt] = extract(obj,x,y)
	nx  = length(x);
	nxc = nx-1;
	if (obj.opt.dischargeisvariable)
		z0  = y(1:nx);
		Q0  = y(nx+1);
		Qt  = reshape(y(nx+2:end),nx,[]);
	else
		y_  = reshape(y,nx,[]);
		z0  = y_(:,1);
		Q0  = obj.Q0_;
		Qt  = y_(:,2:end);
	end
end % extract

