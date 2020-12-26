K=1;

f=@(t,r)oDE(r,K);
[t,Y]=ode45(f,[0 200],0.4);
plot(t,Y)
ylim([0 1])
xlim([0 30])
hold on

function X=oDE(r,K)
X=-(1-K/2)*r-(1/2)*K*r^3;
end

