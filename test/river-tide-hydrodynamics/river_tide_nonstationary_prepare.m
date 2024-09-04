% Tue  6 Sep 14:12:26 CEST 2022
Qlim = [1e3,1e4]
cd   = 2e-3;
W    = 500;
S    = [2.5e-5,5e-5,1e-4];
Q    = formative_discharge(Qlim(1),Qlim(2),'chezy')
h    = normal_flow_depth(Q,W,cd,S,'drag')
L    = h/S;

if (0)
T = (360 + 360)*1440; 
%T = (360 + 20*15)*1440; 
%(20+10)*15*1440;
dt = 60;
r = 0.9907;
%Qs = 0.5*range(Qlim);;
Ts = 360*1440;
%Qmu = Q;
Qmu = mean(Qlim);
Qs  = 0.5*range(Qlim);
%Qs  = 0.5*range(Qmu*[0.2,1.8]);
Qsd = 12*0.25*Q;
p   = 2;
rng(0)
[t,Q,qq] = generate_river_discharge(Qmu,Qsd,Qs,Ts,r,p,T,dt);
plot(t/1440,Q)
t = t/1440;
save('mat/discharge-nonstationary.mat','t','Q');
end
