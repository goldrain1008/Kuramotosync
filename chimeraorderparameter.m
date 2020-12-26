close all
N1=128;
N2=128;
A=0.2;
beta=0.1;
h=0.1;
iter=3000;
omega=0;
theta1=zeros(N1,iter);
theta2=zeros(N2,iter);
r1=zeros(1,iter);
r2=zeros(1,iter);

% IC1=rand(N1,1);
IC1=rand(N1,1);
IC2=rand(N2,1);

% theta1(:,1)=2*asin(sin(IC2));
% theta2(:,1)=2*asin(sin(IC2));
theta1(:,1)=pi*(-1+2*IC1);
theta2(:,1)=pi*(-1+2*IC2);

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
    r1(j)=abs((1/N1)*sum(exp(1i*theta1(:,j))));
    r2(j)=abs((1/N1)*sum(exp(1i*theta2(:,j))));
end

plot(1:iter,r1)
hold on
plot(1:iter,r2)
legend('r1','r2')
ylim([0 1])