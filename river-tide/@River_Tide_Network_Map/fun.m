% Mon 24 Aug 11:01:22 +08 2020
%% morphodynamics of a tidal river
%% either retrive a precomputed scenario or compute and store a new scenario
function [rtn, key, obj] = fun( ...
	          obj ...
		, z1 ...	% [1]     amplitude of incoming wave
		, pz1 ...	% [1]     incoming or total wave
		, omega ...	% [rad/s] anguluar frequency of tide
		, Qmin ...	% [m^3/s] discharge at inflow bc
		, Qmax ...	% [m^3/s] discharge at inflow bc
		, iorder ...
		, S0 ...
		, d_mm ...	% [mm]	  median sediment grain size
		, L ...		% [m]     doamin length
		, pL ...	% [1]     relative location of bifurcation as seen
		     ...	%         from mouth with respect to domain length
		, wa ...	% [m]     width of upstream channel
		, wb ...	% [m]     width of left branch
		, wc ...	% [m]     width of right branch
		, Cd ...	% [1]     drag coefficient
		, ruleQ ...	%  -      rule for discharge division
		, ruleQs ...    %	  rule for sediment division
		, dx ...	% [m]     spatial discretization step
		, Ti ...	% [s]     simulated time span
		, cfl ...	% [1]     cfl condition
		, scheme ... 	%  -      numerical scheme
		)
... % Qs ...	% [kg/s]  sediment transport at inflow bc

%	if (nargin()<15)
%		opt = [];
%	end
	opt = struct();
	opt.xs = 1;

	if (isempty(obj.cmap))
		obj.init();
	end

	key = obj.key(...	
		  z1{1} ...	% [1]     amplitude of incoming wave
		, pz1{1} ...	% [1]     incoming or total wave
		, omega{1} ...	% [rad/s] anguluar frequency of tide
		, Qmin{1} ...	% [m^3/s] discharge at inflow bc
		, Qmax{1} ...	% [m^3/s] discharge at inflow bc
		, iorder{1} ...
		, S0{1} ...
		, d_mm{1} ...		% [mm]	  median sediment grain size
		, L{1} ...		% [m]     doamin length
		, pL{1} ...	% [1]     relative location of bifurcation as seen
		     ...	%         from mouth with respect to domain length
		, wa{1} ...	% [m]     width of upstream channel
		, wb{1} ...	% [m]     width of left branch
		, wc{1} ...	% [m]     width of right branch
		, Cd{1} ...	% [1]     drag coefficient
		, ruleQ{1} ...	%  -      rule for discharge division
		, ruleQs{1} ...    %	  rule for sediment division
		, dx{1} ...	% [m]     spatial discretization step
		, Ti{1} ...		% [s]     simulated time span
		, cfl{1} ...	% [1]     cfl condition
		, scheme{1} ... %  -      numerical scheme
	);
... % Qs ...	% [kg/s]  sediment transport at inflow bc

	% test if simulation was not yet run for the parameter pair
	if (isKey(obj.cmap,key) && ~obj.recompute)
		rtn  = obj.cmap(key);
		t    = rtn.evolution.t;
		zb   = rtn.evolution.zb;
	else
		disp(['recomputing ', key]);

		rtn = River_Tide_Network_2();

		[t, zb] = rtn.evolve_bed_level_scenario(...
		   z1{1}  ...	% [1]     amplitude of incoming wave
		 , pz1{1} ...	% [1]     reflected wave factor
		 , omega{1} ...	% [rad/s] anguluar frequency of tide
		 , Qmin{1} ...	% [m^3/s] discharge at inflow bc
		 , Qmax{1} ...	% [m^3/s] discharge at inflow bc
		 , iorder{1} ... %
		 , S0{1} ...	% [1]	  upstream slope
		 , d_mm{1} ...	% 
		 , L{1} ...	% [m]     doamin length
		 , pL{1} ...	% [1]     relative location of bifurcation as seen
		     ...	%         from mouth with respect to domain length
		 , [wa{1},wb{1},wc{1}] ...	% [m]     width of upstream channel
		 , Cd{1} ...
		 , ruleQ{1} ...	%  -      rule for discharge division
		 , ruleQs{1} ...	%  -      rule for sediment division
		 , dx{1} ...	% [m]     spatial discretization step
		 , Ti{1} ...	% [s]     simulated time span
		 , cfl{1} ...	% [1]     cfl condition
		 , scheme{1} ... 	%  -      numerical scheme
		 , opt ...
	);

		obj.cmap(key) = rtn;
		obj.recflag = true;
		if (obj.autosave)
			obj.save();
		end
	end % if ~iskey
end % River_Tide_Network_Map/fun

