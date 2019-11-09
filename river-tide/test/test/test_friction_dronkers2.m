% seems as if there is something wrong
syms q0 q1 h0 Q0 Q1; f=friction_dronkers((q0 + q1)/h0,q0/h0,Q1/h0,[],true); expand(f)
 syms q0 q1 h0 Q0 Q1 a1; f=friction_trigonometric_dronkers([q0, q1, 0]/h0,[0 a1 0],q0/h0,Q1/h0,[],true); expand(f)

