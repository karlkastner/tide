% Thu  1 Jun 11:29:16 CEST 2017
%% read TPXO output into tidetable object
function [ts, obj] = from_tpxo(obj,filename)
	t = readtable(filename,'delimiter','\t','ReadVariableNames',false);
	%t = readtable('/home/pia/phd/src/tide/mat/tide_besar_mouth/tidal-constituents.dat','delimiter','\t','ReadVariableNames',false);
	l = readtable([ROOTFOLDER,'/src/tide/tidal-constituents-legend.csv'],'delimiter','\t','ReadVariableNames',false);
	n1 = length(l.Var1);
	n2 = length(t.Var2);
	ts = struct();
	ts.tidecon = NaN(max(n1,n2),4);
	ts.freq    = NaN(max(n1,n2),1);
	ts.tidecon(1:n2,1) = t.Var2;
	ts.tidecon(1:n2,3) = t.Var3;
	ts.freq(1:n2) = 1./t.Var1;
	name    = repmat({''},max(n1,n2),1);
	name(1:n1) = l.Var1;
	ts.name = reshape(sprintf('%4s',name{:}),4,[])';
	obj.t = ts;
end

