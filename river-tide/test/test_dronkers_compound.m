% 2017-03-28 16:07:46.956082739 +0200

t=linspace(0,1,10000)';
y = [sin(2*pi*t)  sin(4*pi*t)];
y(:,3)=y(:,1)+1/2*y(:,2);
y=[y(:,1),y(:,3)];
r=range(y);
mr = midrange(y);
% same tidal range
y=bsxfun(@plus,bsxfun(@times,bsxfun(@minus,y,mr),1./r),mr);
u=cdiff(y)/(t(2)-t(1)); plot(u.*abs(u)), mean(u.*abs(u))
