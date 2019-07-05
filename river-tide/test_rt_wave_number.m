% 2018-04-06 18:41:39.385521996 +0200

syms W H g cd omega az1 positive;
Q0 = 0;
Qt=0;
az1 = 0;

[rk, r, k, r_lin] = rt_wave_number(Q0,W,H,cd,omega,az1,Qt,true)
