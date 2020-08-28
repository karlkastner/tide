% Sun  8 Oct 13:08:39 CEST 2017
% Karl Kastner, Berlin
%
%% transform arbitrary to cs-integrated discharge boundary condition
function obj = bc_transformation(obj)

	% for each channel
	for cdx=1:size(obj.bc,3)
	    n_z0 = 0;
	    n_Q0 = 0;			
	
	    % for left and right end of current channel
	    for bdx=1:2
	        % mean water level and discharge
	        switch (obj.bc(bdx,1,cdx).var)
	        case {'z','z0'}
	            n_z0                     = n_z0+1;
	            obj.bc(bdx,1,cdx).set     = true;
	            obj.bc(bdx,1,cdx).var     = 'z';
	        case {'Q','Q0'}
	            %if (~obj.opt.dischargeisvariable)
	                n_Q0         = n_Q0+1;
	                obj.bc(bdx,1,cdx).set = false;
	                obj.bc(bdx,1,cdx).var = 'Q';
	            %else	
	            %end
	        case {''}
	            % nothing to do
	            obj.bc(bdx,1,cdx).set = false;
	        otherwise
	            error('bcfun');
	        end % switch
	        if (obj.bc(bdx,1,cdx).p ~= 1)
	            error('only dirichlet condition supported for mwl so far');
	        end
	        
	        k = 2;
	        df = 1;
	
	        % for each tidal frequency component
	        for jd=k:size(obj.bc,2)
	            switch (obj.bc(bdx,jd,cdx).var)
	            case {'z'}
	                % dQ/dx = -1i*o*z
	                w0 = obj.width(cdx,obj.Xi(cdx,bdx));
	                omega_j              =  (jd-df)*obj.omega;
	                obj.bc(bdx,jd,cdx).rhs        = -1i*omega_j*w0*obj.bc(bdx,jd,cdx).rhs;
	                if ( obj.bc(bdx,jd,cdx).p(2) ~= 0)
	                	error('bc of type dz/dx not yet implemented');
	                end
	                p             = [0, obj.bc(bdx,jd,cdx).p(1), 0];
	                obj.bc(bdx,jd,cdx).p   = p;
	                % change type from level to discharge
	                obj.bc(bdx,jd,cdx).var = 'Q';
	            case {'Q'}
	                % nothing to do, native format
	                % bc(bdx,jd).rhs = bc(bdx,jd).rhs;
	            case {'q'}
	                % Q = w*q
	                obj.bc(bdx,jd,cdx).rhs = obj.bc(bdx,jd,cdx).rhs*W0;
	                % change type from q to Q
	                obj.bc(bdx,jd,cdx).var = 'Q';
	            case {'u'}
	                error('u condition not yet implemented');
	            case {''}
	                % nothing to do
	            otherwise
	                error('bc');
	            end % switch
	            obj.bc(bdx,jd,cdx).set = true;
	        end % for jd, frequency components
	    end % for bdx, left and right end of doamin
	
	%	if (~isreal(obj.z0_downstream) || ~isreal(obj.Q0_))
	%	    error('mean water level and mean discharge must be real');
	%	end
	
	%	if (n_z0 ~= 1)
	%	    error('Water level must be specified on opposit ends of the domain');
	%	end
	
	end % for cdx
end % bc_transform
