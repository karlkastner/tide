% Wed 29 Mar 16:13:27 CEST 2017
% Karl Kastner, Berlin
%% initialize with data
%
% TODO merge with tidal_envelope2
% TODO make class
% TODO max and min vel
% TODO range
function obj = init(obj,t0,dt,H,U)
	if (nargin()<5)
		U = [];
	end

	if (isvector(H))
		H = cvec(H);
		U = cvec(U);
	end

	Ti = obj.Ti;

%	dt   = time(2)-time(1);
	nt   = size(H,1);
	n(1) = round(Ti/dt);
	n(2) = floor(nt/n(1));

	% time at interval centres
	tc = n(1)*dt*((0:n(2)-1)'+0.5);
	
	%dphi_hls = [];
	%t_hls    = [];
	tsh =[];
	tsl =[];
	
	for kdx=1:size(H,2)
		h = H(:,kdx);
	
		% high and low water
		[hw, t_hw, lw, t_lw] = interval_extrema(t0,dt,h,n(1),n(2),obj.order);

		zh(:,kdx)  = hw;
		zl(:,kdx)  = lw;
		tzh(:,kdx) = t_hw;
		tzl(:,kdx) = t_lw;

		if (~isempty(U))
			u = U(:,kdx);
			[uhi, t_uhi, uli, t_uli] = interval_extrema(t0,dt,u,n(1),n(2),obj.order);

			uh(:,kdx)  = uhi;
			ul(:,kdx)  = uli;
			tuh(:,kdx) = t_uhi;
			tul(:,kdx) = t_uli;

			% associated slack water
if (0)
			[t_hws t_lws dt_hws dt_lws] = interval_zeros(t0,dt,u,t_hw,t_lw,n(1),n(2),obj.order);
			
			tsh(:,kdx) = t_hws;
			tsl(:,kdx) = t_lws;
end
			
			%t_hls(:,:,kdx)    = [t_hws, t_lws];
			%dphi_hls(:,:,kdx) = [dphi_hws, dphi_lws];
		end
	end % for each column

	% write back
	obj.zh  = zh;
	obj.zl  = zl;
	obj.tzh = tzh;
	obj.tzl = tzl;
	obj.tsh = tsh;
	obj.tsl = tsl;
	obj.tc  = tc;
	obj.t0 = t0;

	obj.uh  = uh;
	obj.ul  = ul;
	obj.tuh = tuh;
	obj.tul = tul;
	
%	obj.tc=obj.t0 + tc;

	% TODO optionally interpolate to continuous time series
end % function

