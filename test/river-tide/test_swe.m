syms S g cd A w Q dx dw; swe=SWE(); swe.zb=[0; S*dx]; swe.g = g; swe.w=[w;w+dw*dx]; swe.cd = cd; [L V Vi] = swe.fluxmateig([A;Q],w), simplify(V*L*Vi,'ignoreanalyticconstraints',true); [f fA] = swe.flux(0,[A;Q]), swe.source_friction(0,[0;dx],[A;A;Q;Q], [w;w])

