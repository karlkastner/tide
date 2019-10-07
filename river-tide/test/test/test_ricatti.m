function test_ricatti()

	%K2fun = @(x) normpdf(x).^0.1;
	%K2fun = @(x) 1+normpdf(x).^0.5;

	figure(1)
	clf
	I = 2.^(-1:1);
	for idx=1:length(I)

	Xi = idx*[-5,5]*1;
	K2fun = @(x) (1+normcdf(x,0,idx))*(1+I(idx)*i);
	%K2fun  = @(x) 1+exp(-(x-Xi(1))/5);

	%kL    = sqrt(-K2fun(Xi(2)))
	%[x y] = ode45(@odefun,[Xi(2),Xi(1)],kL);
	opt.MaxStep = 1e-1;
	if (1)
		k0    = sqrt(-K2fun(Xi(1)))
		[x y] = ode45(@odefun,[Xi(1),Xi(2)],k0,opt);
	else
	kL    = sqrt(-K2fun(Xi(2)))
	[x y] = ode45(@odefun,[Xi(2),Xi(1)],kL);
	end
	subplot(3,3,3*(idx-1)+1)
	plot(x,abs([(sqrt(-K2fun(x))),y]))
	hold on
	subplot(3,3,3*(idx-1)+2)
	plot(x,imag([(sqrt(-K2fun(x))),y]))
	hold on
	subplot(3,3,3*(idx-1)+3)
	plot(x,real([(sqrt(-K2fun(x))),y]))
	hold on
	end
	
	function kdot = odefun(x,k)
		kdot = -(k.^2 + K2fun(x));
	end
end

