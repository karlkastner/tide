% Sun  2 May 15:49:01 CEST 2021
% Karl Kastner, Berlin
function read_cfg(obj, filename)
	fid = fopen(filename);
	cdx = 0;
	while (1)
		line = fgetl(fid);
		if (~ischar(line))
			break;
		end
		% tokenize
		token = split(line,'=');
		var   = lower(strtrim(token{1}))

		% skip empty lines and comments
		if (~isempty(var) && var(1) ~= '#')

		% reassemble the rest of the line, in case it contains equal signs
		val   = join(token{2:end},'=');
		val   = strtrim(val);

		if (~isempty(val))
			switch (val(1))
			case {'@'}
				% convert string to function handle
				val = eval(val);
			case {'''','"'}
				% leave a string, remove enclosing quotation marks
				val = val(2:end-1);
			case {'0','1','2','3','4','5','6','7','8','+','-',}
				% convert string to number(s)
				val = str2num(val);
			otherwise
				% leave as string
			end
		end % ~isempty(val)

		switch (var)
		case {'omega'}
			obj.rt.omega = val;
		case {'channel'}
			cdx = cdx+1;
			if (1==cdx)
				obj.channel  = River_Tide_Channel();
			else
				obj.channel(cdx)  = River_Tide_Channel();
			end
			obj.channel(cdx).bc = struct();
		case {'length','l'}
			obj.channel(end).L = val;
		case {'nx'}
			obj.channel(end).nx = val;
		% TODO meshing function (similar to set_width)
		case {'width'}
			obj.channel(end).set_width(val);
		case {'zb'}
			obj.channel(end).set_zb(val);
		case {'cd'}
			obj.channel(end).set_cd(val);
		% TODO type of c
		case {'bclvar'}
			obj.channel(end).bc(1).var = val;
		case {'bclval'}
			for idx=1:length(val)
				obj.channel(end).bc(1,idx).rhs = val(idx);
			end
		case {'bclp'} % incoming reflecting (scalar or vector)                         
			obj.channel(end).bc(1).p   = val;
		case {'bclq'} % value vs value of derivative 
			obj.channel(end).bc(1).q   = val;
		case {'bcrvar'}
			obj.channel(end).bc(2,1).var = val;
		case {'bcrval'}
			for idx=1:length(val)
				obj.channel(end).bc(2,idx).rhs = val(idx);
			end
		case {'bcrp'}
			obj.channel(end).bc(2,1).p   = val;
		case {'bcrq'} % value vs value of derivative 
			obj.channel(end).bc(2,1).q   = val;
		otherwise
			error(['Unknown variable ', var]);
		end % switch (var)
		end % if not a comment
	end % until end of file
	neq = size(obj.channel(1).bc,2)
	% copy var, p and q
	for cdx=1:obj.nc
		for bdx=1:2
			for edx=2:neq
				obj.channel(cdx).bc(bdx,edx).var = obj.channel(cdx).bc(bdx,1).var;
				obj.channel(cdx).bc(bdx,edx).p = obj.channel(cdx).bc(bdx,1).p;
				obj.channel(cdx).bc(bdx,edx).q = obj.channel(cdx).bc(bdx,1).q;
			end
			obj.channel(cdx).bc(bdx,edx)
		end
	end
end % read_cfg

