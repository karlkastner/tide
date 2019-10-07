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
%% output :
%% 	p : linear combination of Dirichlet and Neumann boundary condition
%%          1 -> Dirichlet boundary condition
%%	    0 -> Neumann boundary condition
%%      q linear combination of left and right travelling (incoming and outgoing) wave
%% 	rhs = 0 -> homogeneous boundary condition
%%
%%
%% function [rhs, p, q, obj] = bcfun(obj,x,y,ccdx)
function [rhs, p, q, obj] = bcfun(obj,x,y,ccdx)
	Xi = obj.Xi;
	switch (obj.opt.hmode)
	case {'matrix'}
		[rhs, p, q] = bc(x,y,ccdx);
%		switch(ccdx)
%		case {1}
%			%[rhs, p, q] = bc_z0(x,y);	
%		case {2}
%			%[rhs, p, q] = bc_Q1(x,y);
%			[rhs, p, q] = bc(x,y,2);
%		case {3}
%			%[rhs, p, q] = bc_Q2(x,y);
%			[rhs, p, q] = bc(x,y,3);
%		otherwise
%			error('here');
%		end
	otherwise
%		switch (ccdx)
			[rhs, p, q] = bc(x,y,ccdx+1);
%		case {1}
%			[rhs, p, q] = bc_Q1(x,y);
%		case {2}
%			[rhs, p, q] = bc_Q2(x,y);
%		otherwise
%			error('here');
%		end
	end

	function [rhs, p, q] = bc(x,y,id);
		switch (x)
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
	end
%	
%	function [rhs, p, q] = bc_z0(x,y);
%		switch (x)
%		case {Xi(1)}
%			% set value to zero
%			p   = [1,0,0];
%			rhs = 0;
%			q   = [];
%		case {Xi(2)}
%			p = [1,0,0];
%			q = [];
%			% no right bc necessary, this is a first order ode
%			rhs = [];
%		otherwise
%			error('here');
%		end
%	end
%	
%	function [rhs, p, q] = bc_Q1(x,y);
%		switch (x)
%		case {Xi(1)}
%			rhs = obj.bc(1,1).rhs;
%			p   = obj.bc(1,1).p;
%			q   = obj.bc(1,1).q;
%		case {Xi(2)}
%			rhs = obj.bc(2,1).rhs;
%			p   = obj.bc(2,1).p;
%			q   = obj.bc(2,1).q;
%		otherwise
%			error('bcfun');
%		end % switch
%	end
%
%	function [rhs, p, q] = bc_Q2(x,y);
%	switch (x)
%	case {Xi(1)}
%			rhs = obj.bc(1,2).rhs;
%			p   = obj.bc(1,2).p;
%			q   = obj.bc(1,2).q;
%		case {Xi(2)}
%			rhs = obj.bc(2,2).rhs;
%			p   = obj.bc(2,2).p;
%			q   = obj.bc(2,2).q;
%	otherwise
%		error('bcfun');
%	end % switch
%	end
end % River_tide/bcfun

