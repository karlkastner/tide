% Fri Jul 12 17:27:02 UTC 2013
% Karl Kästner, Berlin
%
%% tidal_harmonic analysis
%
%normall null : 
%	all values are taken for one metonic cycle, a period of 19 years
%	LAT  : lowest astronomical tide
%	MLLW : mean astronomical low water tide / mean lower low water (mean of 14 day lows)
%	MLW  : mean low water (mean of semidiurnal lows)
%	NAP  : mean water level in Amsterdam
%	WASG seal level ?
%mean_mw                    = mean(A)
%m =    12*60/15; mean_mlw  = mean(min(reshape(A,n/m, m),[],2))
%m = 14*24*60/15; mean_mllw = mean(min(reshape(A,n/m, m),[],2))
%m                          = min(A)
%lateral flow / transverse / longitudinal / vertical


% t : time in UTC in seconds from year 0
function [A, P, T, legend, Tw, res] = tidal_harmonic_analysis(t, level)

	% load the periods of the tidal constituents, which are constant
	% in space and time
%	T = load([ROOTFOLDER,filesep(),'current/tide/tidal-constituents.csv']);
	addpath([ROOTFOLDER,filesep(),'current/tide/']);
	[T legend Tw] = tidal_constituents();
	T = 3600*T;
	Tw = 3600*Tw;

	% create the regression matrix
	X = zeros(length(t),2*length(T));
	for idx=1:length(T)
		X(:,idx*2-1:idx*2) = [sin(2*pi*t/T(idx)) cos(2*pi*t/T(idx))];
	end

	% regression
	Y = X \ level;

	% TODO improve error estimation
	% cond(X'*X)
	res = X*Y - level;

	% extract amplitude
	A = sqrt(Y(1:2:end).^2 + Y(2:2:end).^2);

	% extract phase
%	P = atan2(Y(1:2:end),Y(2:2:end));
	P = atan2(Y(2:2:end),Y(1:2:end));
%	P = atan(Y(2:2:end)./Y(1:2:end));
%	P = atan(Y(1:2:end)./Y(2:2:end));

end % function tidal_harmonic analysis

