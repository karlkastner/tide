% Tue  8 Oct 11:29:41 PST 2019
% Karl Kästner, Berlin

% river tide series example
% TODO test name = 'friction decomposition';
% TODO test to previous revision

clear
rt_map = River_Tide_Map('./rt-test-map.mat');
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
	, @river_tide_test_12
	};

clear out
nt = length(test_C);
rmse = zeros(nt,2);
fail = false(nt,1);
for idx=1:nt
	testfun = test_C{idx};
	[fail(idx,1),rmse(idx,:),name,out(idx)] = testfun(rt_map,pflag);
	fprintf(['Testing : ', name,'\n']);
	if (fail(idx))
		warning(['Test ',name,' failed']);
	end
end % for idx

rt_map.save();

