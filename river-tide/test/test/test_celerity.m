Q0=[1e3 1e4];
W=500;
H = 16;
cd = 2e-3;
omega=2*pi/86400;
az1=1;
S=16/3e5;
[c c0] = celerity(Q0,W,H,cd,omega,az1,S)

syms Q0 W H cd omega az1 S
[c c0] = celerity(Q0,W,H,cd,omega,az1,S)

expand(simplify(t,'ignoreanalyticconstraints',true))

r = damping_modulus(Q0,W,H,cd,omega,az1); r0=damping_modulus(0,W,H,cd,omega,az1); r=simplify(expand(taylor(r,Q0,0,'order',3)./r0),'ignoreanalyticconstraints',true), subs(r,az1,0)

syms f g; r = simplify(subs(r, Q0, W*sqrt(g*H)*f),'ignoreanalyticconstraints',true); simplify(taylor(r,f,0,'order',3),'ignoreanalyticconstraints',true)

