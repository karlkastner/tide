% Wed 11 Oct 10:18:54 CEST 2017
% TODO alternatively, keep all variables in a single vector and just once compute the indices
function [z0, Q0, Qt] = extract(obj,x,y)
	nx  = length(x);
	nxc = nx-1;
%	switch (obj.opt.hmode)
%	case {'matrix'}
		if (obj.opt.dischargeisvariable)
			z0  = y(1:nx);
			Q0  = y(nx+1);
			Qt  = reshape(y(nx+2:end),nx,[]);
			%z0c = yc(1:nxc);
			%Qtc = reshape(yc(nx+2:end),nxc,[]);
		else
			y_  = reshape(y,nx,[]);
			z0  = y_(:,1);
			Q0  = obj.Q0_;
			Qt  = y_(:,2:end);
			%yc_ = reshape(y,nxc,[]);
			%z0c = yc_(:,1);
			%Qtc = yc_(:,2);
		end
%	otherwise
%		error('here');
%	end % switch
end
