% Fri 13 May 14:44:15 CEST 2016
% Karl Kastner, Berlin
%% import TPXO data into tidetable object 
function tidetable = import_tpxo(tidetable, name)
	dat = load(name);
	if (size(dat,2) < 8)
		error('File empty');	
	end

	tidetable.time  = datenum(dat(:,5), dat(:,3), dat(:,4), dat(:,6), dat(:,7), dat(:,8));
	tidetable.level = dat(:,9);
	tidetable.ux    = 0.01*dat(:,21);
	tidetable.uy    = 0.01*dat(:,22);
	tidetable.umag  = sqrt(tidetable.ux.^2 + tidetable.uy.^2);
	dir             = atan2(tidetable.uy,tidetable.ux);
	dir             = (dir > 0).*dir + (dir < 0).*(dir + 2*pi);
	tidetable.udir  = rad2deg(dir);
	tidetable.dt    = tidetable.time(2) - tidetable.time(1);
end

