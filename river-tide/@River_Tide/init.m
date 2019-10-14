% Sun  8 Oct 13:08:39 CEST 2017
%% provide initial condition by solving the backwater equation for surface level
%% TODO this should not be solved as a ivp but included in the bvp iteration
%% TODO generate the mesh here and precompute fixed values instead of passing functions
%% TODO Q0 should not be a function
%% function obj = init(obj, Xi)
% function obj = init(obj, z_downstream, Q0, Xi)
function obj = init(obj)
	Cfun = @(x) drag2chezy(obj.fun.cd(x));

	obj.bc_transformation();

	obj.fun.Q0        = @(x) obj.Q0_*ones(size(x));

	if (isempty(obj.fun.z0))
		x            = mesh1(obj.Xi,obj.opt.nx,obj.opt.xs);
		obj.backwater.sopt.InitialStep = (x(2)-x(1));
		% TODO no magic numbers
		obj.backwater.sopt.RelTol = 1e-4;
		Qt = 0; % for the initial condition
		[x_, h0_, z0_] = obj.backwater.solve( ...
						obj.Q0_, Qt, Cfun, ...
						obj.fun.width, obj.fun.zb, ...
						obj.z0_downstream(1), obj.Xi);
		obj.tmp.x   = x_;
		obj.tmp.h0  = @(x) interp1(x_,h0_,x,obj.opt.imethod);
		obj.tmp.z0  = @(x) interp1(x_,z0_,x,obj.opt.imethod);
	else
		obj.tmp.x  = linspace(obj.Xi(1),obj.Xi(2),obj.opt.nx).';
		obj.tmp.z0 = obj.fun.z0;
		obj.tmp.h0 = @(x) obj.fun.z0(x) - obj.fun.zb(x);
	end
end % River_Tide/init

