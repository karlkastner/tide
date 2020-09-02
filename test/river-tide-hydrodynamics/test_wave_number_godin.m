% wave number as derived by godin,
% and order of magnitude analysis

K = 2.5e-3
u0 = 1
H = 10
g = 9.81
% for diurnal
omega=2*pi/86400

aa = (3*K*u0^2/H^2)^2
k = sqrt(aa + 8*g*1i*omega*K*u0/H^2)

aa = 0
k = sqrt(aa + 8*g*1i*omega*K*u0/H^2)
