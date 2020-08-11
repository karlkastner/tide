% Tue  8 Oct 11:29:41 PST 2019
% Karl Kästner, Berlin

% river tide test case batch run

% TODO test name = 'friction decomposition';
% TODO test to previous revision
% TODO test case, swap left right                                        
% TODO test case: where does the wobble come from?                              
% TODO scale RMSE2 by max of o2, verify h in backwater computation  

%rt_map = River_Tide_Map('./rt-test-map.mat');
rt_map = River_Tide_Map('/tmp/rt-test-map.mat');
rt_map.recompute = true;
pflag = 1;

test_C = { ... 
	  @river_tide_test_01 ...
	, @river_tide_test_02 ...
	, @river_tide_test_03 ...
	, @river_tide_test_04 ...
	, @river_tide_test_05 ...
	, @river_tide_test_06 ...
	, @river_tide_test_07 ...
	, @river_tide_test_08 ...
	, @river_tide_test_09 ...
	, @river_tide_test_10 ...
	, @river_tide_test_11 ...
	, @river_tide_test_12 ...
	, @river_tide_test_13
	};

nt = length(test_C);
rmse = zeros(nt,2);
fail = false(nt,1);
tab = table();
for idx=1:nt
	disp(['test ', num2str(idx)]);
	tab.Id{idx,1} = idx;
	try
		testfun = test_C{idx};
		[fail(idx,1),rmse(idx,:),name,out_] = testfun(rt_map,pflag);
		tab.Result{idx,1}  = double(fail(idx));	
		tab.RMSE1{idx,1}   = sprintf('%0.4f',rmse(idx,1));
		tab.RMSE2{idx,1}   = sprintf('%0.4f',rmse(idx,2));
		tab.KITER{idx,1}   = rt.out.kiter;
		tab.Name{idx}      = name;
		tab
		out(idx,1:length(out_)) = out_;
		fprintf(['Testing : ', name,'\n']);
	if (fail(idx))
		warning(['Test ',name,' failed']);
	end
	catch e
		tab.Result{idx,1} = true;
		e
	end
end % for idx
tab
%disp(fail,rmse])
rt_map.save();

