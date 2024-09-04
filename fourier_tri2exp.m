%	a cos + b sin = 1/2( a exp ix + a exp-ix + -i b exp ix + i b exp(-ix)
%		      = 1/2( (a - i b) exp ix + (a + i b) exp(-ix)
function cexp = fourier_tri2exp(ctri)
	cexp = zeros(size(ctri));
	cexp(:,1) = ctri(:,1);
	for idx=2:2:size(ctri,2)
		cexp(:,idx) = 0.5*(ctri(:,idx) - 1i*ctri(:,idx+1));
		cexp(:,idx+1) = conj(cexp(:,idx));
	end
end

