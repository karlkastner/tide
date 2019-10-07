% Fri 22 Feb 14:02:42 CET 2019
% Karl Kastner, Berlin
function trn = test_tidal_river_network()
c=1;
switch (c)
case {1} % single channel
	L = 1e5;
	%h = 15;
	zb = -15;
	w = 500;
	%u0 = 1;
	C = 60;
	junction_C = {};
	bc = {1,1,'z',0.0,'z',1.0,0.0;
              1,2,'Q',1e4,'z',0.0,0.0};
case {2} % one channel split in two
	% there is a problem here:
	% 8 unknowns s1l,s1r,c1l,c1r,s2l,slr,c2l,c2r
	% 6 equations:
	% s1(0)
	% c1(0)
	% s2(L2)
	% c2(L2)
	% s1(L1) = s2(0)
	% c1(L1) = c2(0)
	% continuity of the discharge is necessary -> Q1 = Q2
	% continuity of c1 necessary?
	% but why? when con
	% d/dx s1(L1) = d/dx s2(0)
	% d/dx c1(L1) = d/dx s2(0)	

	L = [1e5,1e5];
	h =  15*[1,1];
	w = 500*[1,1];
	u0 = 1*[1,1];
	C = 60;
	bc = [1,1, 1.0,0.0;
              2,2, 0.0,0.0];
	junction = {[1,2;
                     2,1]};
case {3} % 1 channel splitting in two
% we can also introduce the auxiliary variable s_j, c_j -> all branches are equal to sj, cj
% and then sum Q = 0

	s = 0.15;
	L = 10*[s*1e5,2e5,2e5];
	h = [15,15,15];
	w = [500,500/2,500/2];
	u0 = [1,1,1]/8;
	C  = 1e4;
	bc = [1,1,1.0,0;
              2,2,0,0;
	      3,2,0,0];
	junction={[
                   2,1;
                   3,1;
		   1,2;
		]};
case {4} % three channel kapuas like
	% L3 is dummy
	L = [  1e5,   1e5, 2e5];
	h = [   20,    15,  20];
	w = [  500,   220, 500];
	C = 60;
	u0 = [1,1,1]/8;

	junction = { [
		      1,2;
                      2,2;
                      3,1;
		     ]
                    };

	% TODO bc should also contain z0
	bc       = [1,1, 0.0, 1.0;
                    2,1, 0.0, 1.0;
                    3,2, 0.0, 0.0]; 
end

	trn = Tidal_River_Network();
	trn.init(L,zb,w,C,bc,junction_C);
	trn.solve();

	trn.plot_water_level_amplitude();

%	[z,Q] = tidal_river_network()
%	k = [k1, k2, k3];
%	r = [r1, r2, r3];
end

