% Sun 19 Apr 23:43:49 +08 2020
function meta = river_tide_test_metadata()
	% model for river tide
	opt.model_str = 'wave';
	% solver of boundary value problem
%	opt.solver = @bvp2c;

	% number of points along channel
	opt.nx     = 100;

	% change of distance between points along channel 
	opt.xs     = 1; 
	opt.sopt.maxiter = 20;
	
	meta.opt = opt;
end

