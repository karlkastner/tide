% Fri  2 Sep 16:19:44 CEST 2022
function meta = river_tide_nonstationary_metadata()
	meta.mapname_str = 'mat/river-tide-map-nonstationary.mat';
	meta.rtspec_str = 'river-tide-nonstationary.csv';
	meta.folder.d3d  = 'mat/river-tide-nonstationary-';


	% model for river tide
	opt.model_str = 'wave';
	% solver of boundary value problem
%	opt.solver = @bvp2c;

	% number of points along channel
	%opt.nx     = 100;

	% change of distance between points along channel 
	opt.xs     = 1; 
	opt.sopt.maxiter = 200;
	opt.friction_order = 5;
	opt.ode.advective_acceleration = false;
	opt.odefunk = 'odefunk_1';


	
	meta.opt = opt;

	% d3d model parameters
	param = struct();
	param.mdf.Tstop = 2*1440*28;
	% TODO based on cfl and dx
	param.mdf.Dt     = 10;      % minutes
	param.mdf.Roumet = 'C';
	param.mdf.Sub2       = '   ';   % no transport
	param.mdf.Bdf        = 'N'; % no bedform prediction
        param.mor.MorUpd     = false;  % no morupd
        param.mor.IUnderLyr  = 1;     % one well mixed layer
        param.mor.TTLForm    = 1;     % constant transport-layer
	param.mor.ThTrLyr    = 1.0;   % transport-layer thickness
	param.mor.UpdBaseLyr = 4;     % constant base layer
	param.mor.TranspType = 0; % output in kg/s
	%param.mdf.TraFrm = '#vr1993.tra#'
	param.mdf.TraFrm     = 'eh.tra'
	param.mdf.Trtrou     = 'N'
	param.mdf.dt_map     = 30;

	param.bndisharmonic = [true,true];

	meta.param_silent = param;
	meta.param = struct();

end

