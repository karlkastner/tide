t=innerspace(0,1,1e5)';
 Q1=sin(2*pi*t);
 z=sin(2*pi*t);
 p=[];
 for idx=0:5;
 for jdx=1;
 p(idx+1,jdx) = mean(Q1.^idx.*z);
  end;
 end,  p(abs(p)<1e-4) = NaN;
 p*16
sz=1;
 s=1;
 mean((1+s*Q1).^5.*(1-5*z*sz+0*15*z.^2)-0*(1+Q1).^5),
rt_transport_stokes(s,sz,1)

if (0)
	syms e; taylor(1/(1+e)^5,e,0,'order',5), e=(0.1:0.1:0.4)'; o=4:-1:0; s=fliplr([70,-35,15,-5,1].*e.^o), [NaN, fliplr(o); [e, 1./(1+e).^5 - cumsum(s,2)]]
end

