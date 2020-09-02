% Mon  4 Mar 13:59:57 CET 2019

c = 2;
switch (c)
case {1} % single channel
	L = 3e5;
	%h = 15;
	% zb = -15;
	zb = -10;
	w = 500;
	C = 60;
	junction_C = {};

	% channel id, end point id, z/Q, z0, zs1, zc1
	bc = {1,1,'z0',0,;
              1,2,'Q0',1e4};
%              1,2,'z0',10};
case {2}
	% pseudo drawdown discretization
	% TODO, pseudo uniform, pseudo backwater
	nx = 100;
	L_ =  3e5;
	x = L_*linspace(0,1,nx-1)';
	L = L_/(nx-1)*ones(nx-1,1);
	z0_2 = 10;
	c2 = 'uniform';
	c2 = 'backwater';
	switch (c2)
	case {'uniform'}
		zb = -15 + x/L_*z0_2;
	case {'drawdown'}
		zb = -15*ones(nx-1,1);
	case {'backwater'}
		zb = 2*-15 + x/L_*(z0_2+15);
	end
	bc = {   1,1,'z0',0,
	      nx-1,2,'z0',z0_2};
	w = 500; w = w*ones(nx-1,1);
	C = 60; C=C*ones(nx-1,1);
	junction_C = {};
	for idx=2:nx-1
		junction_C{idx-1} =  [idx-1,2 ;
                                     idx,1] ;
	end % for idx
case {3}
	% one channel split in two
	L = 1e5*[1,1,1];
	zb = -15*ones(3,1);
	w = [250,250,500];
	C = 60*ones(3,1);
	bc = {1,1,'z0',0,
	      2,1,'z0',0,
	      3,2,'z0',10};
	junction_C = {[1,2
                       2,2
                       3,1]};
otherwise
	error('here');
end

trn = Tidal_River_Network();
trn.omega = [];
trn.init(L,zb,w,C,bc,junction_C);
trn.solve();
trn.plot_mean_water_level();

