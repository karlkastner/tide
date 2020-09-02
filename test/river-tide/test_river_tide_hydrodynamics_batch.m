% Tue  8 Oct 11:29:41 PST 2019
% Karl Kästner, Berlin

% river tide test case batch run


% TODO mor/network : compare and return test result
% TODO network, mor in batch
% TODO test case : friction decomposition
% TODO test case : comparison to previous revision -> in each test case, load last
% TODO test case : where does the wobble come from?                              
% TODO scale RMSE2 by max of o2, verify h in backwater computation  
% TODO test case mor : run with manning, and compare to seminara
% TODO test case mor : convergence test, compare different schemes
% TODO test case hyd : convergence (z0,z1,Q1)
function [rt_map,rt,tab,out] = test_river_tide_hydrodynamics_batch(pflag,recompute)

	if (nargin()<1)
		pflag = 0;
	end
	if (nargin()<2)
		recompute = true;
	end
	
	% create output file labelled with revision-number
	[retval,rev_str] = system('svn info --show-item revision');
	rev_str = chomp(rev_str);
	filename_str = ['mat/test-rt-hydrodynamics-',rev_str,'.mat'];
	rt_map  = River_Tide_Hydrodynamics_Map(filename_str);
	rt_map.recompute = recompute;
	
	test_C = { ... 
	 	  @test_river_tide_hydrodynamics_01 ...
	 	, @test_river_tide_hydrodynamics_02 ...
	 	, @test_river_tide_hydrodynamics_03 ...
	 	, @test_river_tide_hydrodynamics_04 ...
	 	, @test_river_tide_hydrodynamics_05 ...
	 	, @test_river_tide_hydrodynamics_06 ...
	 	, @test_river_tide_hydrodynamics_07 ...
	 	, @test_river_tide_hydrodynamics_08 ...
		, @test_river_tide_hydrodynamics_09 ...
		, @test_river_tide_hydrodynamics_10 ...
		, @test_river_tide_hydrodynamics_11 ...
		, @test_river_tide_hydrodynamics_12 ...
		, @test_river_tide_hydrodynamics_13
		};
	
	nt = length(test_C);
	clear out;
	clear rt;
	%rmse = zeros(nt,2);
	%fail = false(nt,1);
	tab = table();
	for idx=1:nt
		tab.Id{idx,1} = idx;
		try
			testfun            = test_C{idx};
			[out(idx),rt_]     = testfun(rt_map,pflag);
			disp(['Test #', num2str(out(idx).id), ' ', out(idx).name]);
			
			tab.Result{idx,1}  = double(out(idx).result);	
			tab.RMSE1{idx,1}   = sprintf('%0.4f',out(idx).rmse(1));
			tab.RMSE2{idx,1}   = sprintf('%0.4f',out(idx).rmse(2));
			%tab.KITER{idx,1}   = rt.out.kiter;
			tab.Name{idx}      = out(idx).name;
			rt(idx,1:length(rt_)) = rt_;
			%fprintf(['Testing : ', name,'\n']);
			disp(tab(idx,:));
		if (out(idx).result)
			warning(['Test ',num2str(out(idx).id),' failed']);
		end
		catch e
			tab.Result{idx,1} = true;
			disp(e);
			for sdx=1:length(e.stack)
				disp(e.stack(sdx))
			end
		end
	end % for idx
	disp(tab);
	
	% store result
	rt_map.save();
	% append test result to the map-file
	save(rt_map.filename,'-append','tab','out');
	
end % test_river_tide_hydronynamics_batch

