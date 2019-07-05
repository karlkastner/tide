% Fri 22 Feb 14:41:57 CET 2019
%
%% tide in a fluvial delta channel network, extension of 1D river tide
%% the network is a directed graph
%
classdef Tidal_River_Network < handle
	properties
		% channel length
		L 

		% channel bed level
		zb

		% channel width
		width

		% channel conductivity (Chezy)
		C

		% river flow velocity per channel
		%u0

		% TODO, 0, 1, 2
		omega = 2*pi/86400;

		% boundary conditions
		bc

		% junction conditions
		junction_C

		% channel water level point indices
		channel
	
		% wave number per channel
		k

		% damping rate per channel
		r

		% mean discharge per channel
		Q0
		% mean water level at junction and terminal points
		z0
		% non-zero frequency componets at junctions and terminal points
		z1s
		z1c

		% discretisation matrix and rhs
		A
		rhs
	
		% maximum number of pade iterations
		maxiter = 100;

		g = Constant.gravity;

		% plot distance
		dx=1e3;
	end % properties

	methods
		function obj = Tidal_River_Network()
		end % Tidal River Network
		
		function init(obj,L,zb,width,C,bc,junction_C)
			if (nargin()<4 || isempty(C))
				C = 60*ones(size(L));
			end
			if (isscalar(C))
				C = C*ones(size(L));
			end
			% TODO check for equal length
			obj.L  = L;
			obj.zb = zb;
			obj.width  = width;
			obj.C  = C;
			obj.bc         = bc;
			obj.junction_C = junction_C;

			obj.extract_channels();

			% allocate and initialize
			obj.Q0 = 1; %sqrt(eps)*randn(obj.nc,1);
			obj.z0 = zeros(obj.np,1);

			% TODO allocate higher frequencies

			% TODO test that exactly two conditions are specified per channel
			% TODO test neq and condition number
			% error('end point id must be either 1 (left) or 2 (right)');
			% TODO check range of bdx
		end % init

		function extract_channels(obj)
			channel = zeros(obj.nc,2);

			% generate channel matrix
			for idx=1:obj.nb
				% water level points that are end points
				id = obj.bc{idx,1};
				jd = obj.bc{idx,2};
				channel(id,jd) = idx;
			end % idx

			% for each junction
			% water level points that are junctions
			for jdx=1:obj.nj
				junction = obj.junction_C{jdx};
				% for each connected channel
				for cdx=1:size(junction,1);
					% channel id
					cid = junction(cdx,1);
					% end point id
					eid = junction(cdx,2);
					if (channel(cid,eid) == 0)
						channel(cid,eid) = obj.nb+jdx;
					else
						error('channel with duplicate end points');
					end % if
				end % for cdx
			end % for jdx

			obj.channel = channel;
		end % extract_channel

		% pseudo members
		% number of channels (graph edges)
		function nc = nc(obj)
			nc = length(obj.L);
		end % n

		% number of terminal points
		function nb = nb(obj)
			nb = size(obj.bc,1);
		end

		% number of junctions
		function nj = nj(obj)
			nj = length(obj.junction_C);
		end

		% number of water level points (graph vertices)
		function np = np(obj)
			np = obj.nb + obj.nj;
		end
	
		function nf = nf(obj)
			nf = length(obj.omega);
		end

		function cd = cd(obj)
			cd     = sqrt(obj.g)./obj.C;
		end

		function h = channel_depth(obj)
			z0 = obj.z0;
			zb = obj.zb;
			channel = obj.channel;
			h  = [z0(channel(:,1)) - zb, ...
			      z0(channel(:,2)) - zb];
		end % channel_depth
	end % methods
end % Tidal_River_Network

