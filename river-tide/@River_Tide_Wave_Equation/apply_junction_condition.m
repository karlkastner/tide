% Wed 23 Feb 11:56:47 CET 2022
% TODO, this is for the explicit solver, implement for implicit solver
function Q = apply_junction_condition(obj,t,Q)
	for jdx=1:length(obj.junction)
		% channel index
		cid = obj.junction(jdx).cid;
		% endpoint index
		eid = obj.junction(jdx).eid;

		% global index of junction grid point
		% (note that each connecting channel has an own grid point at the junction)
		at = [];
		% global index of first grid poin in the channel from the junction
		in = [];

		% for each connected channel
		for cdx=1:length(cid)
			% global index of junction point amd index of last but junction point
			if (1 == eid(cdx))
				% junction is at first end of channel
				at(cdx) = obj.first(cid(cdx));
				in(cdx) = at(cdx)+1;
			else
				% junctio is at second end of channel
				at(cdx) = obj.last(cid(cdx));
				in(cdx) = at(cdx)-1;
			end
		end % for idx

		% width of channels at junctions
		% TODO the width should be a function and evaluated at end points
		w      = obj.w(cid);

		% extract distance from grid point to point at junction from first grid point in the channel
		dx = obj.x(at) - obj.x(in);

		% extract value at first grid point in channel
		Qin = Q(in);
	
		% set up linear system
		A   = zeros(length(cid));
		rhs = zeros(length(cid),1); 
		% first row : sum(Q) = 0
		A(1,:) = 1e-6; % this is better conditioned than 1
		% remaining rows:
		% TODO this only sets dz/dt not z!
		% z1 = z2 -> 1/w dQ/dx|_channel 1 = 1/w dQ/dx|_channel 2
		% this is ill conditioned for Qt=0 (as dQ_river/dt ~ 0)
		for idx=1:length(cid)-1
			A(idx+1,idx)   = +1/w(idx)*1/dx(idx);
			A(idx+1,idx+1) = -1/w(idx+1)*1/dx(idx+1);
			rhs(idx+1)     = -1/w(idx)*1/dx(idx)*Qin(idx) + 1/w(idx+1)*1/dx(idx+1)*Qin(idx+1); 
		end % for idx
w
format long
A
rhs
		% solve for discharge at junction
		Qat = A \ rhs;
rcond(A)
rank(A)
condest(A)
Qin
Qat
dx
A.*Qat.'
pause
		% write to global indices
		Q(at) = Qat;
	end % for jdx (each junction)
end % apply_junction_condition

