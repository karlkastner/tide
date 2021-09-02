% Wed  9 Oct 15:23:10 PST 2019
function rt = hydrodynamic_scenario(rt_map,zl,ql,qr,zb,Q0,w0,Cd,omega,Lx,opt,Tseason)

	opt.oflag = [true(1,length(zl)-1),false(1,4-length(zl))];

	bc          = struct();

	% for each frequency component
	for idx=1:length(zl)
		% boundary at left end
		% mean sea level
		bc(1,idx).var = {'z'};
		bc(1,idx).rhs = zl(idx);
		% Dirichlet condition
		if (idx==1)
			bc(1,idx).p   = 1;
		else
			bc(1,idx).p   = [1,0];		
		end
		% wave entering from left
		bc(1,idx).q   = [1,ql];
	
		% boundary condition at right end
		bc(2,idx).var = {'z'};
		bc(2,idx).rhs = 0;
		% Dirichlet condition
		if (idx==1)
		bc(2,idx).p   = 1;
		else
		bc(2,idx).p   = [1,0];
		end
		bc(2,idx).q   = [qr,1];
	end

	for idx=length(zl)+1:4
		for jdx=1:2
		    bc(jdx,idx).var = {};
		    bc(jdx,idx).rhs = [];
		    bc(jdx,idx).p   = [];
		    bc(jdx,idx).q   = [];
		end
	end

	% mean discharge
	% (overwrite right end bc for zero-frequency component)
	% river discharge

	bc(2,1).var = {'Q'};
	if (nargin()>11 && ~isempty(Tseason))
		bc(2,1).Qseason = Q0;
		bc(2,1).Tseason = 60*Tseason;
		bc(2,1).phase_season = 0;
	else
		bc(2,1).rhs = Q0;
	end
	% domain size
	Xi        = [0, Lx];

	% solve with model
	% TODO, this does not account for seasonal variation (!)
	rt  = rt_map.fun({Xi} ... % Q0,
			, {w0}, {Cd}, {zb}, omega ...
			,  bc(1,1).var, {bc(1,1).rhs}, {bc(1,1).p} ...
			,  bc(2,1).var, {bc(2,1).rhs}, {bc(2,1).p} ...
			,  bc(1,2).var, {bc(1,2).rhs}, {bc(1,2).p}, {bc(1,2).q} ...
			,  bc(2,2).var, {bc(2,2).rhs}, {bc(2,2).p}, {bc(2,2).q} ...
			,  bc(1,3).var, {bc(1,3).rhs}, {bc(1,3).p}, {bc(1,3).q} ...
			,  bc(2,3).var, {bc(2,3).rhs}, {bc(2,3).p}, {bc(2,3).q} ...
			, opt);
	rt.channel(1).bc = bc;
end

