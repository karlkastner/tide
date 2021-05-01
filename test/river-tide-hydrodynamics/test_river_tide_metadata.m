% Sun 19 Apr 23:43:49 +08 2020
function meta = river_tide_test_metadata()
	[retval,rev_str] = system('svn info --show-item revision');
	rev_str      = chomp(rev_str);
	meta.mapname_str = ['mat/test-rt-hydrodynamics-',rev_str,'.mat'];
	meta.folder.d3d  = 'mat/test-river-tide-';
	meta.testspec    = [ROOTFOLDER,'/src/lib/tide/test/river-tide-hydrodynamics/test-river-tide.csv'];

	% model for river tide
	opt.model_str = 'wave';
	% solver of boundary value problem
%	opt.solver = @bvp2c;

	% number of points along channel
	opt.nx     = 100;

	% change of distance between points along channel 
	opt.xs     = 1; 
	opt.sopt.maxiter = 20;
	opt.sopt.maxiter = 200;
	opt.friction_order = 5;
	opt.ode.advective_acceleration = true;
	opt.ode.advective_acceleration = false;
	opt.odefunk = 'odefunk_1';
	
	meta.opt = opt;

	% d3d model parameters
	param = struct();
	param.mdf.Tstop = 2*1440*28;
	% TODO based on cfl and dx
	param.mdf.Dt     = 10;      % minutes
	param.mdf.Roumet = 'C';
	param.mdf.Sub2   = '   ';   % no transport
	param.mdf.Bdf        = 'N'; % no bedform prediction
        param.mor.MorUpd     = false;  % no morupd
        param.mor.IUnderLyr  = 1;     % one well mixed layer
        param.mor.TTLForm    = 1;     % constant transport-layer
	param.mor.ThTrLyr    = 1.0;   % transport-layer thickness
	param.mor.UpdBaseLyr = 4;     % constant base layer
	param.mor.TranspType = 0; % output in kg/s
	%param.mdf.TraFrm = '#vr1993.tra#'
	param.mdf.TraFrm = 'eh.tra'
	param.mdf.Trtrou = 'N'
	param.mdf.dt_map = 60;
	param.bndisharmonic = true;

	meta.param_silent = param;
	meta.param = struct();
end

