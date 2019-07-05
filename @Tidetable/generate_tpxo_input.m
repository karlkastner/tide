% Sun May 19 12:12:49 UTC 2013
% todo, lea years

%% generate tpxo input table
%% Note: superseeded by perl script

function A = generate_table()
	lat = 0;
	long = 109;
	n = 1;
	%       j  f  m  a  m  j  j  a  s  o  n  d
	dom = [31 28 31 30 31 30 31 31 30 31 30 31];
%	mih = 0;
%	mih=[0 15 30 45];
	mih = [0 5 10 15 20 25 30 35 40 45 50 55];

	A = zeros(365*24*12,8);

	% year
	for y=2013:2014
	% month
	for m=1:12
	% day
	for d=1:dom(m)
	% hour
	for h=0:23
	% minutes
		for mi=mih;
			A(n,:) = [ lat long y m d h mi 0];
			n = n+1;
		end
	end % h
	end % d
	end % m
	end % y
	A = A(1:n-1,:);
	f=fopen('pontianak/lat_lon_time','w');
	fprintf(f,'%d %d %d %d %d %d %d %d\n',A');
	fclose(f);
end

