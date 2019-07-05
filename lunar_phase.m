% Sat 12 May 15:13:46 CEST 2018
%% lunar phase
% c.f. wikipedia
function [nmod,new,full,n] = lunar_phase(t)
	% dt = 29.530588853; t0 = datenum('01/01/1900');
	
	% first new month in current era
	c0 = 5.597661 + datenum('2000-01-01 00:00:00');
	% length synodic month
	c1 = 29.5305888610;
	% drift of length of the synodic month
	c2 = 1.02026e-10;
	% t = c0 + c1*n + c2*n^2
	t = cvec(t);
	o = ones(size(t));
	n = roots2(fliplr([c0-t,c1*o,c2*o]));
	n = n(:,1);
	% fraction 
	nmod = mod(n,1);
	% mark new
	n_   = [nmod(1)-1.0; nmod];
	%new  = (n_(1:end-1) < 0.0) & (n_(2:end) >= 0.0);
	new  = (n_(1:end-1) > 0.5) & (n_(2:end) <  0.5);
	full = (n_(1:end-1) < 0.5) & (n_(2:end) >= 0.5);
end

