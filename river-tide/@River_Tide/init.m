% Sun  8 Oct 13:08:39 CEST 2017
%% solve backwater equation for surface level
%% TODO this should not be solved as a ivp but included in the bvp iteration
%% TODO generate the mesh here and precompute fixed values instead of passing functions
function obj = init(obj, z_downstream, Q0, Xi)

	C = @(x) sqrt(obj.g/obj.cdfun(x));
	
	obj.Xi   = Xi;
	obj.Q0fun  = @(x) Q0; %*ones(size(x));
	obj.z_downstream = z_downstream;

	if (isempty(obj.fun.z0))
		% TODO apply xs
		%x = obj.mesh1d(
		%obj.backwater.sopt.InitialStep = (Xi(2)-Xi(1))/obj.opt.nx;
		x            = mesh1(Xi,obj.opt.nx,obj.opt.xs);
		obj.backwater.sopt.InitialStep = (x(2)-x(1));
		obj.backwater.sopt.RelTol = 1e-4;
		[x_, h0_, z0_] = obj.backwater.solve(Q0,0,C,obj.wfun,obj.fun.zb,obj.z_downstream(1),obj.Xi);
		obj.tmp.x   = x_;
		obj.tmp.h0  = @(x) interp1(x_,h0_,x,obj.opt.imethod);
		obj.tmp.z0  = @(x) interp1(x_,z0_,x,obj.opt.imethod);
	else
		obj.tmp.x  = linspace(obj.Xi(1),obj.Xi(2),obj.opt.nx).';
		obj.tmp.z0 = obj.fun.z0;
		obj.tmp.h0 = @(x) obj.fun.z0(x) - obj.fun.zb(x);
	end
end % RT::init

