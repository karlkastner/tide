% Mon  4 Nov 19:15:55 +08 2019
%	(a + bi)*exp(it)  = a cos x + b i cos x  + a i sin x + b i^2 sin x
%	(a - bi)*exp(-it) = a cos x - b i cos x  - a i sin x + b i^2 sin x
%			    2 a cos x - 2 b sin x
function ctri = fourier_exp2tri(c_exp)
	ctri = zeros(size(c_exp));
	ctri(:,1) = c_exp(:,1);
	for idx=2:2:size(c_exp,2)
		ctri(:,idx)   =  2*real(c_exp(:,idx));
		ctri(:,idx+1) = -2*imag(c_exp(:,idx));
	end
end

