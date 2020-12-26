close all
N1=180;
N2=180;
A=0.2;
beta=0.1;
h=0.1;
iter=3000;
% pd=makedist('tLocationScale','mu',0.5,'sigma',0.0007,'nu',1);
% omega=random(pd,N1,1);
omega=0.5;
theta1=zeros(N1,iter);
theta2=zeros(N2,iter);
theta1(:,1)=zeros(N1,1);
theta2(:,1)=ICchimera;
% theta1(:,1)=zeros(N1,1);
% theta2(:,1)=asin(sin(randn(N2,1)))+0.2079;
% 2*asin(sin(randn(N2,1)));
% IC1=rand(N1,1);
% IC2=rand(N2,1);
% theta1=zeros(N1,iter);
% theta2=zeros(N2,iter);
% theta1(:,1)=pi*(-1+2*IC1);
% theta2(:,1)=pi*(-1+2*IC2);

for j=1:iter
    k1=h*chimerapop1(theta1(:,j),theta2(:,j),A,beta,N1,N2,omega);
    l1=h*chimerapop2(theta1(:,j),theta2(:,j),A,beta,N1,N2,omega);
    k2=h*chimerapop1(theta1(:,j)+0.5*k1,theta2(:,j)+0.5*l1,A,beta,N1,N2,omega);
    l2=h*chimerapop2(theta1(:,j)+0.5*k1,theta2(:,j)+0.5*l1,A,beta,N1,N2,omega);
    k3=h*chimerapop1(theta1(:,j)+0.5*k2,theta2(:,j)+0.5*l2,A,beta,N1,N2,omega);
    l3=h*chimerapop2(theta1(:,j)+0.5*k2,theta2(:,j)+0.5*l2,A,beta,N1,N2,omega);
    k4=h*chimerapop1(theta1(:,j)+k3,theta2(:,j)+l3,A,beta,N1,N2,omega);
    l4=h*chimerapop2(theta1(:,j)+k3,theta2(:,j)+l3,A,beta,N1,N2,omega);
    theta1(:,j+1)= theta1(:,j)+(1/6)*(k1+2*k2+2*k3+k4);
    theta2(:,j+1)= theta2(:,j)+(1/6)*(l1+2*l2+2*l3+l4);
    
%     x=cos(theta1(:,j));
%     y=sin(theta1(:,j));
%     s=linspace(-pi,pi,100);
%     cx=cos(s);
%     cy=sin(s);
%     plot(x,y,'ro',cx,cy);
%     axis([-1 1 -1 1]);
%     axis square
%     xlabel('Population 1')
%     title(['t= ' num2str(j)])
%     
%     subplot(1,2,2)
%     u=cos(theta2(:,j));
%     v=sin(theta2(:,j));
%     plot(u,v,'kx',cx,cy);
%     axis([-1 1 -1 1]);
%     axis square
%     xlabel('Population 2')
%     title(['t= ' num2str(j)])
%     

end
% plotting order parameter
Z1=(1/N1)*sum(exp(1i*theta1));
Z2=(1/N2)*sum(exp(1i*theta2));
r1=abs(Z1);
r2=abs(Z2);
phi1=angle(Z1);
phi2=angle(Z2);
psi=phi1-phi2;
% plot(r2)
% hold on
% plot(r1)
% xlim([0 iter])
% ylim([0 1])
% title(['A= ' num2str(A)])

figure(15)
plot(r2)
xlim([0 iter])
ylim([0 1])
xlabel('time')
ylabel('r')
title('A=0.2')
% plot(r2.*cos(psi),r2.*sin(psi))
% xlim([0 1])
% ylim([-0.5 0.5])
% hold on
% plot(r2(1)*cos(psi(1)),r2(1)*sin(psi(1)),'rx')

