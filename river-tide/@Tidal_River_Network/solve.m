% Fri 22 Feb 14:01:09 CET 2019
%
%% solve for the tide in a fluvial chanel network
%%
%% boundary condition at end points not connected to junctions
%% 	[ channel 1 id, endpoint id (1 or 2), s0, c0
%%        ...
%%         channel n id, endpoint id (1 or 2), s0, c0]
%%
%% conditions at junctions are specified as cells
%%	 each cell contains an nx2 array
%%	 n : number of connecting channels
%%	 [channel id1, endpoint id (1 or 2), ...
%%         channel idn, endpoint id (1 or 2)]
%%
%% every tidal species for each channel has 4 unknowns
%% these are 2x2 unknowns for the sin + cos of left and right going wave
%
% TODO 0	switch to complex analysis 
% TODO 10	friction depending on tidal and river flow
% TODO 15	feed back of tide on mwl
% TODO 20	change to mid channel (-dx,+dx)
% TODO 40	make matrix sparse (setup through buffer) 
% TODO 40	merge with River_Tide (in particular, the mwl computation not any more through ode23)
% TODO 50	second species and overtides
% TODO 75	tidal flats
% TODO 100	convergence (w,h)
%
function [z,Q] = solve(obj)

	%n = obj.n;

	% allocate memory
	if (0 == obj.nf)
		ne = obj.np + obj.nc;
	else
		ne = obj.np + obj.nc*(1 + 4*obj.nf);
	end
	A   =  zeros(ne);
	rhs = ne;

	% TODO pade iteration
	sol_ = 0;
	iter = 0;
	while (true)
		if (obj.nf > 0)
			[r, k] = damping_modulus_river(obj.Q0,obj.width,obj.channel_depth,obj.cd,obj.omega);

			obj.r = -r;
			obj.k = -k;
		end

		[A,rhs] = matrix0(obj,A,rhs);

		if (obj.nf > 0)
			[A,rhs] = matrix1(obj,A,rhs);
		end

		obj.A   = A;
		obj.rhs = rhs;
		
		% solve
		sol = A \ rhs;

		p    = 0.5;
		sol  = p*sol + (1-p)*sol_;
		sol_ = sol;
		
		% extract result and store
		% inside loop, because this is required for the matrix setup

		% water level per junction and end point
		obj.z0 = sol(1:obj.np);

		% discharge per channel
		obj.Q0 = sol(obj.np+1:obj.np+obj.nc);

%		size(A)
%		[A,rhs]
%		obj.z0
%		obj.Q0
%		pause
	

		if (obj.nf > 0)
			% store
			nc = obj.nc;
			np = obj.np;
			sol = sol(nc+np+1:end);
			obj.z1s = reshape(sol(1:2*nc),2,[])';
			obj.z1c = reshape(sol(2*nc+1:end),2,[])'; 
		end
		
		iter = iter+1;
		if (iter >= obj.maxiter)
			break;
		end
	end

end % solve

% matrix main tidal species
function [A,rhs] = matrix1(obj,A,rhs)
	% allocate memory
	% this is a sparse matrix, but as it is small, we do not care
%	A   = zeros(4*n);
%	rhs = zeros(4*n,1);

	nc = obj.nc;
	
	bc         = obj.bc;
	junction_C = obj.junction_C;

	% TODO, these have to be adapted iteratively
	k = obj.k;
	r = obj.r;

	L = cvec(obj.L);
	w = obj.width;
	x = [zeros(nc,1),L];

	% skip equations for z0 and Q0
	neq0 = obj.np+obj.nc;
	
	% boundary condition for each outlet/inlet
% TODO change from Q to z
	for bdx=1:obj.nb
		id = bc{bdx,1};
		xi = x(id,bc{bdx,2});
		% zc = -1/(iow) dQc/dx = -1/(iow) e( -r c + k s) hat Qc

		switch (bc{bdx,5})
		case {'z'}
			% sin left going wave
			A(neq0+2*bdx-1,neq0+2*id-1)     = exp(-r(id)*xi)*cos( k(id)*xi); 
			A(neq0+2*bdx-1,neq0+2*id-1+2*nc) = exp(-r(id)*xi)*sin(+k(id)*xi);

			% sin right going wave
			A(neq0+2*bdx-1,neq0+2*id)       = exp(-r(id)*(L(id)-xi))*cos( k(id)*(L(id)-xi));
			A(neq0+2*bdx-1,neq0+2*id+2*nc)   = exp(-r(id)*(L(id)-xi))*sin(+k(id)*(L(id)-xi));

			% cos left going
			A(neq0+2*bdx,neq0+2*id-1)       = exp(-r(id)*xi)*sin(-k(id)*xi);
			A(neq0+2*bdx,neq0+2*id-1+2*nc)   = exp(-r(id)*xi)*cos( k(id)*xi);

			% cos right going
			A(neq0+2*bdx,neq0+2*id)         = exp(-r(id)*(L(id)-xi))*sin(-k(id)*(L(id)-xi));
			A(neq0+2*bdx,neq0+2*id+2*nc)     = exp(-r(id)*(L(id)-xi))*cos(+k(id)*(L(id)-xi));
		case {'Q'}
			% sin left going wave
			A(neq0+2*bdx-1,neq0+2*id-1)     = -1/(1i*o*w(id))*exp(-r(id)*xi) ...
						          *( -r(id)*cos(k(id)*xi) + k(id)*sin(k(id)*xi) );
			A(neq0+2*bdx-1,neq0+2*id-1+2*nc) = -1/(1i*o*w(id))*exp(-r(id)*xi) ...
							*(-r(id)*sin(+k(id)*xi) + k(id)*cos(+k(id)*xi));

			% sin right going wave
			A(neq0+2*bdx-1,neq0+2*id)       = -1/(1i*o*w(id))*exp(-r(id)*(L(id)-xi)) ...
								* (r(id)*cos(+k(id)*(L(id)-xi)) + k(id)*sin(k(id)*(L(id)-xi)));
			A(neq0+2*bdx-1,neq0+2*id+2*nc)   = -1/(1i*o*w(id))*exp(-r(id)*(L(id)-xi)) ...
								* (r(id)*sin(+k(id)*(L(id)-xi)) - k(id)*cos(k(id)*(L(id)-xi)));

			% cos left going
			A(neq0+2*bdx,neq0+2*id-1)       = -1/(1i*o*w(id))*exp(-r(id)*xi) ...
							  *(-r(id)*sin(-k(id)*xi) - k(id)*cos(-k(id)*xi));
			A(neq0+2*bdx,neq0+2*id-1+2*nc)   = -1/(1i*o*w(id))*exp(-r(id)*xi) ...
							  *(-r(id)*cos( k(id)*xi) - k(id)*sin(k(id)*xi));
			
			% cos right going
			A(neq0+2*bdx,neq0+2*id)   	= -1/(1i*o*w(id))*exp(-r(id)*(L(id)-xi)) ... 
							  *(r(id)*sin(-k(id)*(L(id)-xi)) -k(id)*cos(-k(id)*(L(id)-xi)));
			A(neq0+2*bdx,neq0+2*id+2*nc)     = -1/(1i*o*w(id))*exp(-r(id)*(L(id)-xi)) ...
							  *(r(id)*cos(+k(id)*(L(id)-xi)) + k(id)*sin(+k(id)*(L(id)-xi)));
		otherwise
			error('here');
		end % switch boundary condition type

		% right hand side for sin
		rhs(neq0+2*bdx-1) = bc{bdx,6}; 

		% right hand side for cos
		rhs(neq0+2*bdx) = bc{bdx,7};
	end

	% junctions 2*(n-1)*n conditions
	neq=neq0+2*obj.nb;
	for jdx=1:obj.nj
		junction = junction_C{jdx};
		nj       = size(junction,1);

		% identity of the level z1lr_i = z1lr_1
		c1 = junction(1,1);
		x1  = x(c1,junction(c1,2))
		for idx=2:nj
			c2 = junction(idx,1)
			x2 = x(c2,junction(idx,2));

			% sin(omega t)
			neq=neq+1;

			% channel 1 : left going
			A(neq,neq0+2*c1-1)           =  1/(1i*o*w(c1))*exp(-r(c1)*x1) ...
							*( -r(c1)*cos(+k(c1)*x1) - k(c1)*sin(+k(c1)*x1));
			A(neq,neq0+2*c1-1+2*nc)      =  1/(1i*o*w(c1))*exp(-r(c1)*x1) ... 
							*( -r(c1)*sin( k(c1)*x1) + k(c1)*cos(+k(c1)*x1) );

			% channel 1 : right going
			A(neq,neq0+2*c1)             =   1/(1i*o*w(c1))*exp(-r(c1)*(L(c1)-x1)) ...
							*(+r(c1)*cos(+k(c1)*(L(c1)-x1)) + k(c1)*cos(+k(c1)*(L(c1)-x1)));
			A(neq,neq0+2*c1+2*nc)        =   1/(1i*o*w(c1))*exp(-r(c1)*(L(c1)-x1)) ...
							*(+r(c1)*sin(+k(c1)*(L(c1)-x1)) -k(c1)*cos(+k(c1)*(L(c1)-x1)));

			% channel 2 : left going
			A(neq,neq0+2*c2-1)           =  -1/(1i*o*w(c2))*exp(-r(c2)*x2) ...
							*( -r(c2)*cos(+k(c2)*x2) - k(c2)*sin(+k(c2)*x2));
			A(neq,neq0+2*c2-1+2*nc)      =  -1/(1i*o*w(c2))*exp(-r(c2)*x2) ... 
							*( -r(c2)*sin( k(c2)*x2) + k(c2)*cos(+k(c2)*x2) );

			% channel 2 : right going
			A(neq,neq0+2*c2)             =   -1/(1i*o*w(c2))*exp(-r(c2)*(L(c2)-x2)) ...
							*(+r(c2)*cos(+k(c2)*(L(c2)-x2)) + k(c2)*cos(+k(c2)*(L(c2)-x2)));
			A(neq,neq0+2*c2+2*nc)        =   -1/(1i*o*w(c2))*exp(-r(c2)*(L(c2)-x2)) ...
							*(+r(c2)*sin(+k(c2)*(L(c2)-x2)) -k(c2)*cos(+k(c2)*(L(c2)-x2)));

			% right hand side is zero

			% cos(omega t)
			neq=neq+1;

			% channel 1 : left going
			A(neq,neq0+2*c1-1)      =  1/(1i*o*w(c1))*exp(-r(c1)*x1) ...
						   *(-r(c1)*sin(-k(c1)*x1) -k(c1)*cos(-k(c1)*x1));
			A(neq,neq0+2*c1-1+2*nc) =  1/(1i*o*w(c1))*exp(-r(c1)*x1) ...
						   *(-r(c1)*cos(+k(c1)*x1) - k(c1)*sin(+k(c1)*x1));

			% channel 1 : right going
			A(neq,neq0+2*c1)        =   1/(1i*o*w(c1))*exp(-r(c1)*(L(c1)-x1)) ...
						    *(+r(c1)*sin(-k(c1)*(L(c1)-x1)) + cos(-k(c1)*(L(c1)-x(1))));
			A(neq,neq0+2*c1+2*nc)   =   1/(1i*o*w(c1))*exp(-r(c1)*(L(c1)-x1)) ...
						    *(+r(c1)*cos(+k(c1)*(L(c1)-x1)) + k(c1)*sin(+k(c1)*(L(c1)-x1)));

			% channel 2 : left going
			A(neq,neq0+2*c2-1)      =  -1/(1i*o*w(c2))*exp(-r(c2)*x2) ...
						   *(-r(c2)*sin(-k(c2)*x2) -k(c2)*cos(-k(c2)*x2));
			A(neq,neq0+2*c2-1+2*nc) =  -1/(1i*o*w(c2))*exp(-r(c2)*x2) ...
						   *(-r(c2)*cos(+k(c2)*x2) - k(c2)*sin(+k(c2)*x2));

			% channel 2 : right going
			A(neq,neq0+2*c2)        =   -1/(1i*o*w(c2))*exp(-r(c2)*(L(c2)-x2)) ...
						    *(+r(c2)*sin(-k(c2)*(L(c2)-x2)) + cos(-k(c2)*(L(c2)-x(1))));
			A(neq,neq0+2*c2+2*nc)   =   -1/(1i*o*w(c2))*exp(-r(c2)*(L(c2)-x2)) ...
						    *(+r(c2)*cos(+k(c2)*(L(c2)-x2)) + k(c2)*sin(+k(c2)*(L(c2)-x2)));

			% right hand side is zero
		 end % for idx

		% discharge condition
		neq = neq+1;
		for idx=1:nj
			ci  = junction(idx,1);
			xi  = x(ci,junction(idx,2));
			
			% sin omega t, sum sQ1l_i + sQ1r_i == 0
			% left going
			A(neq,neq0+2*ci-1)        =  exp(-r(ci)*xi)*cos(+k(ci)*xi);
			A(neq,neq0+2*ci-1+2*nc)   =  exp(-r(ci)*xi)*sin(+k(ci)*xi);
			% right going
			A(neq,neq0+2*ci)          =  exp(-r(ci)*(L(ci)-xi))*cos(+k(ci)*(L(ci)-xi))
			A(neq,neq0+2*ci  +2*nc)   =  exp(-r(ci)*(L(ci)-xi))*sin(+k(ci)*(L(ci)-xi));
			
			% sum cQ1l_i + cQ1r_i == 0
			% left going
			A(neq+1,neq0+2*ci-1)      =  exp(-r(ci)*xi)*sin(-k(ci)*xi);
			A(neq,neq0+2*ci-1+2*nc)   =  exp(-r(ci)*xi)*cos(+k(ci)*xi);
			% right going
			A(neq+1,neq0+2*ci)        =  exp(-r(ci)*(L(ci)-xi))*sin(-k(ci)*(L(ci)-xi));
			A(neq+1,neq0+2*ci+2*nc)   =  exp(-r(ci)*(L(ci)-xi))*cos(+k(ci)*(L(ci)-xi));
		end
		neq = neq+1;
		% rhs is zero
	end % for jdx
end % matrix1

% matrix for mwl and mean flow
function [A, rhs] = matrix0(obj,A,rhs)
	% number of water level points
	np = obj.np;

	% number of channels
	nc = obj.nc;

	% allocate memory
	%A   = zeros(np+nc);
	%rhs = zeros(np+nc,1);

	zb  = obj.zb;
	z0  = obj.z0;
	Q0  = obj.Q0;
	cd  = obj.cd;
	w   = obj.width;
	%x   = obj.x;
	dx  = obj.L; %x(:,2)-x(:,1);
	h   = obj.channel_depth();
	hc  = mean(h,2);
	% wide channel approximation
	Rh = hc;
	% TODO, influence of tide
	friction_term = cd.*abs(Q0)./(Rh.*w.^2.*hc.^2);

	% for each channel, set stationary flow equation (chezy)
	% g dz/dx + cd Q0^2/(w^2*h^3)
	ne = 0;
	for cdx=1:nc
		ne       = ne+1;
		id       = obj.channel(cdx,1);
		jd       = obj.channel(cdx,2);
		A(ne,id) = -obj.g/dx(cdx);
		A(ne,jd) = +obj.g/dx(cdx);
		% mid point rule
		A(ne,np+cdx) = -friction_term(cdx);
	end 
	
	% junctions
	for jdx=1:obj.nj
		junction = obj.junction_C{jdx};
		nj      = size(junction,1); 
		p1      = obj.channel(junction(1,1),junction(1,2));
if (0)
		% equality of water level in all channels
		% p1 == pi, automatically quaranteed
		for idx=2:nj
			ne = ne+1;
			pi = obj.channel(junction(idx,1),junction(idx,2));
%p1
%pi
%pause
			% water level : z_j = z_1
			A(ne,p1) =  1;
			A(ne,pi) = -1;
			% rhs(ne) is zero
		end % for idx
end

		% discharge at junctions adds to zero
		ne = ne+1;
		for idx=1:nj
			cid=junction(idx,1);
			% in or outflow (
			switch (junction(idx,2))
			case {1}
				A(ne,np+cid)  =  -1;
			case {2}
				A(ne,np+cid) =  +1;
			otherwise
				error('here');
			end
			% rhs(ne) is zero
		end % for idx
	end % for idx
	
	% boundary condition
	for idx=1:obj.nb
		bc  = obj.bc(idx,:);
		cid = bc{1};
		pid = obj.channel(cid,bc{2});
		ne  = ne+1;
		switch (bc{3})
			case {'z'}
				A(ne,pid)  = 1;
				rhs(ne,1)  = bc{4};
			case {'Q'}
				A(ne,np+cid) = 1;
				rhs(ne,1)    = bc{4};
			otherwise
				error('here');
		end % switch
	end % for idx
	%pause
end % matrix0

