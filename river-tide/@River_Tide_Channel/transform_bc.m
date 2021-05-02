% Sun  8 Oct 13:08:39 CEST 2017
% Karl Kastner, Berlin
%
%% transform arbitrary to cs-integrated discharge boundary condition
function transform_bc(obj)

	omega = obj.rt.omega;

	% for each channel
	    n_z0 = 0;
	    n_Q0 = 0;			
	
	    % for left and right end of current channel
	    for bdx=1:2
	        % mean water level and discharge
	        switch (obj.bc(bdx,1).var)
	        case {'z','z0'}
	            n_z0                     = n_z0+1;
	            obj.bc(bdx,1).set     = true;
	            obj.bc(bdx,1).var     = 'z';
	        case {'Q','Q0'}
	            %if (~obj.opt.dischargeisvariable)
	                n_Q0         = n_Q0+1;
	                obj.bc(bdx,1).set = false;
	                obj.bc(bdx,1).var = 'Q';
	            %else	
	            %end
	        case {''}
	            % nothing to do
	            obj.bc(bdx,1).set = false;
	        otherwise
	            error('bcfun');
	        end % switch
	        if (obj.bc(bdx,1).p ~= 1)
	            error('only dirichlet condition supported for mwl so far');
	        end
	        
	        k = 2;
	        df = 1;

		Xi = [0,obj.L];	
	        % for each tidal frequency component
	        for jd=k:size(obj.bc,2)
	            switch (obj.bc(bdx,jd).var)
	            case {'z'}
	                % dQ/dx = -1i*o*z
	                w0 = obj.width(Xi(bdx)); %obj.Xi(bdx));
	                omega_j              =  (jd-df)*omega;
	                obj.bc(bdx,jd).rhs        = -1i*omega_j*w0*obj.bc(bdx,jd).rhs;
	                if ( obj.bc(bdx,jd).p(2) ~= 0)
	                	error('bc of type dz/dx not yet implemented');
	                end
	                p             = [0, obj.bc(bdx,jd).p(1), 0];
	                obj.bc(bdx,jd).p   = p;
	                % change type from level to discharge
	                obj.bc(bdx,jd).var = 'Q';
	            case {'Q'}
	                % nothing to do, native format
	                % bc(bdx,jd).rhs = bc(bdx,jd).rhs;
	            case {'q'}
	                % Q = w*q
	                obj.bc(bdx,jd).rhs = obj.bc(bdx,jd).rhs*W0;
	                % change type from q to Q
	                obj.bc(bdx,jd).var = 'Q';
	            case {'u'}
	                error('u condition not yet implemented');
	            case {''}
	                % nothing to do
	            otherwise
	                error('bc');
	            end % switch
	            obj.bc(bdx,jd).set = true;
	        end % for jd, frequency components
	    end % for bdx, left and right end of doamin
	
	%	if (~isreal(obj.z0_downstream) || ~isreal(obj.Q0_))
	%	    error('mean water level and mean discharge must be real');
	%	end
	
	%	if (n_z0 ~= 1)
	%	    error('Water level must be specified on opposit ends of the domain');
	%	end
	
end % transform_bc

