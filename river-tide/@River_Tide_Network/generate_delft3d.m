% Tue  5 Nov 20:05:15 +08 2019
% Karl Kastner, Berlin
%
%% generate a Delft3D 4 model for the channel network
%
% note : amplitude should be specified as sin (purely imaginary z10), so that
%        model starts up with a flat surface
%
% note : tidal amplitude to depth ration in frictionless cases should not exceed
%	 1%, as otherwise the nonlinear steepening of the wave by advective
%	 acceleration causes very high spurious oscillations in Delft3D,
%	 as the numerical sheme is not TVD
%
function d3d = generate_delft3d(obj, folder, param, param_silent, opt);
	param_silent.omega = obj.rt.omega;

	if (obj.nc>1)
		error('multiple channels are not yet supported')
	end

	if (~isfieldorprop(opt,'Lc'))
		opt.Lc = 0;
	end
	if (~isfieldorprop(opt,'Cz_end'))
		opt.Cz_end = 50;
	end

	nn  = [2, obj.hydrosolver.nx];
	Xi  = obj.hydrosolver.xi;

	mesh = StructuredMesh();
	mesh.generate_rectangle(  0.5*[1,-1] ...
			            , Xi + sqrt(eps) ...
				    , [nn(1),nn(2)]);

	% along-channel coordiantes for gradually varying meshes
	x          = obj.hydrosolver.out(1).x;
	% avoid zero, as delft3d interprets it as NaN
	x(x==0)    = sqrt(eps);
	% extend domain
	dx = x(end)-x(end-1);
	if (opt.Lc-x(end)>=dx)
		x  = [rvec(x),x(end)+(dx:dx:(opt.Lc-x(end)))];
		
	else
		opt.Cz_end = [];
	end
	%y     = repmat(rvec(x),2,1);
	y = [1/2;-1/2]*ones(1,length(x));

	% set width, X is here the accross and Y the along channel coordinate
	w0         = rvec(obj.channel(1).width(cvec(x)));
	%X          = mesh.X;
%	x          = x;
	Y          = w0.*y;
	%d3d.mesh.Y = Y;

	% roughness
	%Cd  = obj.cd(1,x);
	%Cz  = drag2chezy(Cd);
	%Cz  = min(1e4,Cz);
	if (~isempty(opt.Cz_end))
		% TODO no magic numbers
		Czmin = 1e3;
		x_  = [cvec(obj.channel(1).x); opt.Lc]
		Cz_ = [cvec(min(Czmin,drag2chezy(obj.channel(1).cd))); opt.Cz_end];
		param_silent.Cz = @(x,y) interp1(x_, Cz_, y,'linear','extrap');
	else
		param_silent.Cz = @(x,y) interp1(obj.channel(1).x,min(drag2chezy(obj.channel(1).cd),1e3),y,'linear','extrap');
	end

	% must follow computation of Cz
	mesh.X     = Y;
	mesh.Y     = repmat(rvec(x),2,1);

	% TODO quick fix for off by 1/2 error
%	x  = Y(1,:);
%	zb = obj.zb(1,cvec(x));
%	zb = [zb(2:end); 2*zb(end)-zb(end-1)];
%	zb  = repmat(rvec(zb),size(Y,1),1);
	param_silent.zb = @(x,y) interp1(obj.channel(1).x,obj.channel(1).zb,y,'linear','extrap');

	% boundary condition
	switch (obj.channel(1).bc(1,2).var)
	case {'z'}
		zs0(1) = obj.channel(1).bc(1,1).rhs;
		for idx=2:size(obj.channel(1).bc,2)
			zs0(idx) = 2*obj.channel(1).bc(1,idx).rhs;
		end
	case {'Q'}
		% iow z + dQ/dx = 0
		zs0 = 0;
		for idx=2:size(obj.channel(1).bc,2)
			%zs0 = obj.channel(1).bc(1,idx,1).rhs;
			dQ_dx    =  obj.channel(1).bc(1,idx).rhs;
			zs0(idx) = -dQ_dx./(1i*(idx-1)*obj.rt.omega*w0(1));
		end
	otherwise
		error('here');
	end
	param_silent.zs0 = zs0;

	% this fails when the length of the domain exceeds 5000km
	% due to a segfault in d3d

	% extend the domain at the right end
	% TODO this is inappropriate for test cases with flow from left to right

	%	Lc = opt.Lc;
	% limit chezy coefficients for frictionless cases
	% TODO no magic numbers
%	Cz1 = min(Cz(end),100);
%	Cz1 = min(Cz(end),100);
%		Cz_ = fun(x_,Cz(end),Cz_end);
%		x_  = x_ + x(end);
		% extend domain
%		x    = [ x, x_];
%		Cz   = [Cz,Cz_];
%	else
%		x_ = [];
%		Cz_ = [];
%	end
%	end
%	Y = repmat(x,size(Y,1),1);
%	X = repmat(X(:,1),1,length(x));
	
%	Cz2 = repmat(Cz,size(Y,1),1);
	% TODO, inner2outer
%	Cz2 = [Cz2, Cz2(:,end)];
%	Cz2 = [Cz2; Cz2(end,:)];

	param_silent.z00 = obj.channel(1).bc(1,1,1).rhs;
	param_silent.Q0 = obj.channel(1).bc(2,1,1).rhs;
%	param.nn = nn;

	param_silent.sediment = obj.sediment;

	mesh.generate_delft3d(folder,param,param_silent);

	% smoothly vary between initial and final value
	function y = fun(x,a,b)
		p1 = 0;
		p2 = 1;
		x = (x-x(1))/(x(end)-x(1));
		y(x<=p1) = a;
		fdx = x>p1 & x<p2;
		y(fdx) = exp(log(b) + (log(a)-log(b))*(cos(0.5*pi*(x(fdx)-p1)/(p2-p1))));
		%y(fdx) = ((b) + ((a)-(b))*(cos(0.5*pi*(x(fdx)-p1)/(p2-p1))));
		%y(fdx) = ((b) + ((a)-(b))*0.5*(1+cos(pi*(x(fdx)-p1)/(p2-p1))));
		y(x>=p2) = b;
	end
end % generate_delft3d


