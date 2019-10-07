% Thu 12 Oct 13:29:42 CEST 2017
%% key for storing a scenario
%%
%% function [key obj] = key(obj,varargin)
function [key obj] = key(obj,varargin)
		key = '';
		for idx=1:length(varargin)
			if (ischar(varargin{idx}))
				key = [key, varargin{idx}]
			else 
				if (~isscalar(varargin{idx}))
					error('arguments have to be strings or scalaras');
				end % if
				key = [key, sprintf('%e',varargin{idx})];
			end
			if (idx < length(varargin))
				key = [key,' '];
			end
		end % for idx
end % River_Tide_Map/key

