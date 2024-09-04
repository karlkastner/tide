% Fri 30 Oct 17:11:18 +08 2020
% fourider derivatve
% y     = sum_j hat y_j exp(i j omega t + i k_j x)
% dy/dt = sum_j i j omega exp(i j omega t + i k_j x)
function Q_t = fourier_derivative(obj, Q)
	n     = size(Q,2)-1;
	k     = 0:n; 
	omega_k = rvec(obj.omega(k));
	Q_t   = (1i*omega_k).*Q;
end

