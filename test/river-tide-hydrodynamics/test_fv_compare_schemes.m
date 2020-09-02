% 2017-11-14 23:45:17.885372831 +0100

% https://depts.washington.edu/clawpack/clawpack-4.3/book/chap13/swhump1/

% hump and dambreak
% double break

% TODO odd and even number of steps

L  = 10;
nx = 64;
%128;
h  = 0.5;
a  = 4/100;
g  = 9.81;
Ti = [0 0.3*L/sqrt(g*(0.5*a+h))];

w  = 3;
%x = linspace(-5,5);
ic_hump       = @(x) [h+a*normpdf(w*x); zeros(size(x))];
ic_dambreak   = @(x) 1+1*(x<0);
ic_dambreak2  = @(x) 1+1*(x<-1/6*L)+1*(x<1/6*L);
ic_wavepacked = @(x) 1 + normpdf(0.5*w.*x).*cos(2*pi*x.*w*0.5);
% TODO make wavepacked travelling (subtract left velocity)

% ic = @(x) 1+2*(x<0);
%ic = @(x) 1+abs(x<0.5);
% xtick(-5:5); ylim(0.5*[-1 1]+1)
%ic = ic_wavepacked;
ic = ic_hump;

%x = linspace(-L/2,L/2,1e3);
%plot(x,ic(x))
%pause

nlim = 6;
k = 2;

%  lax friedrich
%- lw + updwind
%- lw + beam-warming
%- lw + minmod
%- lw + superbee
% TODO rename LW into RAE

	limiter_C = {'lax_friedrich',
		     'upwind',
                     'lax_wendroff',
		     'vanLeer',
                     'minmod',
                     'superbee',
                     'monotized_central',
		     'beam_warming',
		     'fromm'
			};

% TODO check convergence (with and w/o source term)

	close all
	Y = [];
	Vol = [];
	figure(1)
	clf
	for idx=1:nlim
%length(limiter_C)
		switch (limiter_C{idx})
		case {'lax_friedrich'}
			fv = Lax_Friedrich();
		otherwise
			fv         = Reconstruct_Average_Evolve();
			%FL = Flux_Limiter();
			fv.limiter = @(varargin) Flux_Limiter.(limiter_C{idx})(varargin{:});
		end
		fv.bcfun{1} = @SWE.bc_nonreflecting;
		fv.bcfun{2} = @SWE.bc_nonreflecting;
		fv.icfun    = ic;
		fv.pde = SWE();
		fv.init([-L/2,L/2],nx);
		[T Y_] = fv.solve(Ti);
		Y(:,idx) = Y_(1:nx,end);
		%if (1==idx)
		figure(1)
		plot(interp1(T,Y_(1:nx,:)',linspace(Ti(1),Ti(2),10))');
%		end
		Y_ = interp1(T,Y_(1:nx,:)',linspace(Ti(1),Ti(2),100))';
		Vol(:,idx) = sum(Y_(1:nx,:));
		figure(101)
		subplot(2,3,idx)
		cla
		plot([fv.x,-fv.x],[Y(:,idx),Y(:,idx)]);
	end
	figure(100);
	subplot(2,2,1)
	plot(fv.x,Y);
	legend(limiter_C{1:nlim})
	subplot(2,2,2)
	plot(Vol)
	

if (1)

end
