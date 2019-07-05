% Fri 15 Jul 17:25:26 CEST 2016
% Karl Kastner, Berlin
%% run TPXO to generate time series
function obj = generate(obj,foldername,timeshift,latlon,zone,y0,ye);
	tpxoname   = [foldername,filesep(),'tide_.out'];

	% get path
	w = what('Tidetable');
	% ROOTFOLDER,'/src/lib/tide/@Tidetable/
	p = w.path;

	% generate the TPXO input and call TPXO
	command = [p,filesep(),'tidetable_generate.sh ', foldername ...
				,' ',num2str(latlon(1)),' ',num2str(latlon(2)) ...
				,' ',num2str(y0), ' ', num2str(ye)];
	disp(command);
	[ret, str] = system(command);
	if (ret ~= 0)
		error(['Bash: ',str]);
	else
		str
	end
	
	% read TPXO output
	obj.import_tpxo(tpxoname)

	% load coordinates
%	latlon = load(latlonname);

	% set name
	obj.placename = basename(foldername);
	
	% shift time zone
	obj.time = obj.time + timeshift;

	% set coordinates
	if (nargin > 4)
		[obj.x obj.y] = latlon2utm(latlon(1),latlon(2),zone);
	end
end

