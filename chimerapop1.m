function dXdt=chimerapop1(X,Y,A,beta,N1,N2,Omega)
 mu = (1+A)/2;
 nu = (1-A)/2;
 alpha = (pi/2)-beta;
 dXdt = Omega+(mu/N1)*sum(sin(ones(N1,1)*X'-(X*ones(1,N1))-alpha),2)+(nu/N2)*sum(sin(ones(N1,1)*Y'-(X*ones(1,N2))-alpha),2);
end

