% Wed  9 Oct 15:23:10 PST 2019

function rt = hydrodynamic_scenario(rt_map,zl,ql,qr,zb,Q0,w0,Cd,omega,Lx,opt,Tseason,Qseason,nx,description)

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

	% fill remaining frequency components with zeros
	if (0)
	for idx=length(zl)+1:4
		for jdx=1:2
		    bc(jdx,idx).var = {};
		    bc(jdx,idx).rhs = [];
		    bc(jdx,idx).p   = [];
		    bc(jdx,idx).q   = [];
		end
	end
	end

	% mean (river) discharge
	% (overwrite right end bc for zero-frequency component)
	% TODO make this variable left and right channel
	bc(2,1).var = {'Q'};
	bc(2,1).rhs = Q0;
	omega_season = NaN;

	% domain size
	Xi        = [0, Lx];

	opt.nx = nx;

	% solve with model
	rt  = rt_map.fun({Xi} ... % Q0,
			, {w0}, {Cd}, {zb}, omega ...
			,  bc(1,1).var, {bc(1,1).rhs}, {bc(1,1).p} ...
			,  bc(2,1).var, {bc(2,1).rhs}, {bc(2,1).p} ...
			,  bc(1,2).var, {bc(1,2).rhs}, {bc(1,2).p}, {bc(1,2).q} ...
			,  bc(2,2).var, {bc(2,2).rhs}, {bc(2,2).p}, {bc(2,2).q} ...
			,  bc(1,3).var, {bc(1,3).rhs}, {bc(1,3).p}, {bc(1,3).q} ...
			,  bc(2,3).var, {bc(2,3).rhs}, {bc(2,3).p}, {bc(2,3).q} ...
			, opt);
	% TODO, add to fun as argument
	rt.description = description;

	% the seasonal discharge is not accounted for, but has to be written here,
	% furthermore, the untransformed boundary conditions have to be used
	% to generate the analoguous delft3d-model
	rt.channel(1).bc = bc;
	if (nargin()>11 && ~isempty(Tseason))
		rt.channel(1).bc(2,1).Qseason = Qseason;
		rt.channel(1).bc(2,1).Tseason = Physics.SECONDS_PER_DAY*Tseason;
		rt.channel(1).bc(2,1).phase_season = 0;
	end
end

