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
%	 acceleration causes very high spurious oscillations in Delft3D
%
function d3d = generate_delft3d(obj, tid, param_silent, Lc);

	d3d        = Delft3D();

	% TODO no magic names
	d3d.folder = ['mat/test-river-tide-',num2str(tid)];

	% note that time step in d3d is given in minutes (!)
	% TODO no magic numbers
	dt_bc   = 1/(24);  % days
	itdate  = '2020/01/01';
	tratype = 103;

	% model properties
	z00 = obj.bc(1,1,1).rhs;
	Q0  = obj.bc(2,1,1).rhs;
	Xi  = obj.hydrosolver.xi;
	nn  = [2, obj.hydrosolver.nx];

	mkdir(d3d.folder);
	d3d.read_all(d3d.templatefolder());
	d3d.mdf.set_filenames(d3d.base);
	% d3d.tratype = param.tratype;
	d3d.itdate  = datenum(itdate);
	d3d.tratype = tratype;

	d3d.mesh.generate_rectangle(  0.5*[1,-1] ...
			            , Xi + sqrt(eps) ...
				    , [nn(1),nn(2)]);

	% extend domain by fading out elements (2e6m)
	Y  = d3d.mesh.Y;
	X  = d3d.mesh.X;
	x = Y(1,:);
	dx = x(end)-x(end-1);
	Cd  = obj.cd(1,x);
	Cz  = drag2chezy(Cd);
	Cz  = min(1e4,Cz);

	% this fails when the lenght of the domain exceeds 5000km
	% due to a segfault in d3d
	if (0)
		L = 1e6;
		s = 2;
		xend = x(end);
		x_   = cumsum([logspace(log10(dx),log10(L),round(s*log2(L/dx)))]);
	else
		if (Cz(end) > 100)
			Cz_end = 1;
		else
			Cz_end = Cz(end);
		end
		%L   = 4e6;
		%b   = 1.00125;
		%n   = round(log(1-(1-b)*L/dx)/log(b))
		%x_  = x(end) + cumsum(dx*(b.^(1:n)));
		%Cz_ = interp1([x_(1),x_(end)],[Cz(1),1],x_,'linear');
		%Cz_ = exp(interp1([x_(1),x_(end)],log([Cz(1),Cz_end]),x_,'linear'));
		if (Lc-x(end)>=dx)
		dx = x(2)-x(1);
		x_ = dx:dx:(Lc-x(end));
		%xend = 0;
		%for idx=1:6
		%x_  = [x_, xend+x(2:end)];
		%xend=x_(end);
		%end
		%x_  = [x_, x_(end)+x(2:end)];
		%x_  = [x_, x_(end)+x(2:end)];
		% taper off
		%Cz_ = (Cz(end)-Cz_end)*0.5*(1+cos(2*pi*(x_-x_(1))/(x_(end)-x_(1))))+Cz_end;
		%Cz_(x_>0.5*x_(end)) = Cz_end;
		%Cz_ = exp((log(Cz(end))-log(Cz_end))*0.5*(1+cos(pi*(x_-x_(1))/(x_(end)-x_(1))))+log(Cz_end));
		Cz_ = fun(x_,Cz(end),Cz_end);
		x_ =x_ + x(end);
		else
			x_ = [];
			Cz_ = [];
		end
	end
function y = fun(x,a,b)
	p1 = 0; %0.5;
	p2 = 1;
	x = (x-x(1))/(x(end)-x(1));
	y(x<=p1) = a;
	fdx = x>p1 & x<p2;
	y(fdx) = exp(log(b) + (log(a)-log(b))*(cos(0.5*pi*(x(fdx)-p1)/(p2-p1))));
	%y(fdx) = ((b) + ((a)-(b))*(cos(0.5*pi*(x(fdx)-p1)/(p2-p1))));
	%y(fdx) = ((b) + ((a)-(b))*0.5*(1+cos(pi*(x(fdx)-p1)/(p2-p1))));
	y(x>=p2) = b;
end

	x    = [ x, x_];
	Cz   = [Cz,Cz_];

	Y = repmat(x,size(Y,1),1);
	X = repmat(X(:,1),1,length(x));

	% expand 
%	Cz  = 1e3;
	Cz2 = repmat(Cz,size(Y,1),1);
%	Cz2(Y>xend) = 50;
	Cz2 = [Cz2, Cz2(:,end)];
	Cz2 = [Cz2; Cz2(end,:)];
	d3d.mesh.export_delft3d_rgh([d3d.folder,'/','delft3d.rgh'],Cz2);

	% set width, X is here the accross and Y the along channel coordinate
	w0 = rvec(obj.width(1,cvec(x)));
	X  = w0.*X;
	d3d.mesh.X = X;
	d3d.mesh.Y = Y;

	% TODO quick fix for off by 1/2 error
%	zs = zs';
%	zs = [zs(1,:);
%	      mid(zs);];
%	zs = [mid(zs);
%             zs(end,:)];
%	zs = zs';
	x  = Y(1,:);
	%zb = obj.zb(1,mid(cvec(x)));
	zb = obj.zb(1,cvec(x));
	zb = [zb(2:end); 2*zb(end)-zb(end-1)];
%zb(1)
%2*zb(1)-zb(2)
%pause
	%zb = [(2*zb(1)-zb(2));zb];

	% set bed level
	zb  = repmat(rvec(zb),size(Y,1),1);
	%zb  = reshape(obj.zb(1,Y(:)),size(Y))';
	h0  = -zb(1);

	switch (obj.bc(1,2,1).var)
	case {'z'}
		z10 = obj.bc(1,2,1).rhs;
	case {'Q'}
		% iow z + dQ/dx = 0
		dQ_dx = obj.bc(1,2,1).rhs;
		z10 = -dQ_dx./(1i*obj.omega*w0(1))
	otherwise
		error('here');
	end
	%z10 = 0.1*z10*exp(0.5i*pi);
	d3d.mesh.Z = zb;

	% TODO, set automatically after generation
	nn                      = d3d.mesh.n;
	d3d.mdf.mdf.dat.MNKmax  = sprintf(' %d %d %d',[nn(2)+1,nn(1)+1,1]);

	% initial condition
	zs = z00*ones(size(zb));
	if (abs(Q0)>0)
		zs = max(zs,zb+h0);
	end
%	zs = [mid(zs(1:end-),zs(end)];
	%zs = [2*zs(:,1)-zs(:,2),zs(:,1:end-1)];
	zs = [2*zs(:,1)-zs(:,2),zs]; %;(:,1:end-1)];
	zs = mid(zs')';
	%zs = [zs,2*zs(:,end)-zs(:,end-1)];
	
	% TODO, take filename from mdf
	d3d.write_ini([d3d.folder,'/delft3d.ini'],zs); 

	% configure mdf
	d3d.mdf.set(param_silent.mdf);
	d3d.mor.set(param_silent.mor);
%	d3d.mdf.set(param.mdf);
	% crashes otherwise
	%d3d.mdf.mdf.dat.Filsed = '##'

	bct = struct();

	% boundary condition
	t   = (0:dt_bc:param_silent.mdf.Tstop/1440)';
	z   = z00 + real(z10*exp(1i*obj.omega*t*86400));
	dx  = Y(1,end)-Y(1,end-1)
	S0  = (zb(1,end)-zb(1,end-1))./dx
	w0  = w0(end)

	g = Physics.gravity;
	c = sqrt(g*h0);
	% ignore flow vel, for time being
	% c*dt < dx
	dt_max = dx/c;
	fprintf('dt_max %f\n',dt_max);

% TODO merge bct and bnd -> consistency has to be checked when reading of files
% TODO assign id automatically by index
% "name" is a redundant field
% "type" can be determined outomatically

	% outflow boundary (downstream)
	bct(1).id       = 1;
	bct(1).type     = 'Waterlevel';
	bct(1).location = 'Outflow';
	bct(1).dt_d     = dt_bc;
	bct(1).time     = t;
	bct(1).val      = z;

	% inflow boundary (upstream)
	bct(2).id       = 2;
	bct(2).type     = 'Discharge';
	bct(2).location = 'Inflow';
	bct(2).dt_d     = dt_bc; %1/24;
	bct(2).time     = t; % ^[0,1/24,param.mdf.Tstop/1440];
 			     % Q0*[1e-3,1,1]; 
	S0_  = S0;
	Cd_  = Cd(end);
	Q0_t = discharge_step_response(t*86400,h0,w0,S0_,Cd_);
	bct(2).val      = Q0_t;

	d3d.bct = bct;

	% flow at open boundaries
	bnd          = struct();
	bnd(1).name  = 'Outflow';
	bnd(1).left  = [1,2];
	bnd(1).right = [1,nn(1)];
	bnd(1).type  = 'Z';

	bnd(2).name  = 'Inflow';
	bnd(2).left  = [nn(2)+1,2];
	bnd(2).right = [nn(2)+1,nn(1)];
	bnd(2).type  = 'T';
	d3d.bnd = bnd;

	%gsd_ = d3d.set_fractions(1e3*meta.d50_m,1.0);
	d3d.set_fractions(obj.sediment.d_mm,1.0);

	d3d.default_bcc();

	d3d.write_all();
end % generate_delft3d

