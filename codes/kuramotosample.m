clc, clf, clear, close all
N=180;
K=4;
h=0.1;
iter=1000;
theta=zeros(N,iter);
theta(:,1)=2*pi*rand(N,1);
Omega=randn(N,1);

for j=1:iter
    k1=kuramoto(theta(:,j),K,N,Omega);
    k2=kuramoto(theta(:,j)+0.5*h*k1,K,N,Omega);
    k3=kuramoto(theta(:,j)+0.5*h*k2,K,N,Omega);
    k4=kuramoto(theta(:,j)+h*k3,K,N,Omega);
    theta(:,j+1)= theta(:,j)+(h/6)*(k1+2*k2+2*k3+k4);
    x=cos(theta(:,j));
    y=sin(theta(:,j));
    s=linspace(0,2*pi,100);
    cx=cos(s);
    cy=sin(s);
    plot(x,y,'o',cx,cy);
    axis([-1 1 -1 1]);
    axis square
    title(['t= ' num2str(j)])
    drawnow
end

    
