% Thu  1 Jun 09:12:43 CEST 2017
%% wrapper for TPXO generated tidal time series
classdef T_Tide < handle
	properties
		t
		% reverse index for names
		id
		%freq
		%amplitude
		%serr
		%phase
	end

	methods 
		function obj = T_Tide(t)
			if (nargin() > 0)
			for idx=1:length(t)
				obj(idx).t = t(idx);
				%obj(idx)=struct();
				%obj(idx).amplitude = t(idx).tidecon(:,1);
				% obj(idx).serr.a    = t(idx).tidecon(:,2);
				% obj(idx).phase     = t(idx).tidecon(:,3);
				%obj(idx).serr.p    = t(idx).tidecon(:,4);
				%name               = t(idx).name;
			end % for idx
			end % if
		end % T_Tide

		% pseudo members
		function f = freq(obj)
			f = obj.t.freq;
		end
		function a = amplitude(obj)
			a = obj.t.tidecon(:,1);
		end
		function p = phase(obj)
			p = obj.t.tidecon(:,3);
		end
		function se_a = serr_amplitude(obj)
			se_a = obj.t.tidecon(:,2);
		end
		function se_a = serr_phase(obj)
			se_a = obj.t.tidecon(:,4);
		end
		function name = name(obj)
			name  = obj.t.name();
			name  = regexprep(mat2cell(name,ones(length(name),1),size(name,2)),' *','');
		end

		function obj = sort(obj,field,varargin)
			for idx=1:length(obj)
				% get order of field
				[void sdx] = sort(obj(idx).(field),varargin{:});
				% sort the data fields
				obj(idx).t.tidecon = obj(idx).t.tidecon(sdx,:);
				obj(idx).t.freq    = obj(idx).t.freq(sdx,:);
				obj(idx).t.name    = obj(idx).t.name(sdx,:);
			end % for idx
		end % sort
	end % methods
end % classdef

