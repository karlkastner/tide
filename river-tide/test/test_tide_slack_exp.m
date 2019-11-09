% Tue  5 Nov 09:55:58 +08 2019
n = 6;
c = rand(n,3) + 1i*[zeros(n,1),rand(n,2)];

[r0,y0] = tide_slack_exp(c);
t       = linspace(0,1);

[r_,y_] = tide_low_high_exp(c);

clf
for idx=1:n
subplot(2,3,idx)
y       = c(idx,1) + c(idx,2)*exp(2i*pi*t) + c(idx,3)*exp(4i*pi*t);
y       = real(y);
plot(t,y)
hold on
plot(r0(idx,:),y0(idx,:),'o');
hline(0)
plot(r_(idx,:),y_(idx,:),'*');

end
