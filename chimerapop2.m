function dYdt=chimerapop2(X,Y,A,beta,N1,N2,Omega)
 mu = (1+A)/2;
 nu = (1-A)/2;
 alpha = (pi/2)-beta;
 dYdt = Omega+(mu/N2)*sum(sin(ones(N2,1)*Y'-Y*ones(1,N2)-alpha),2)+(nu/N1)*sum(sin(ones(N2,1)*X'-Y*ones(1,N1)-alpha),2);
end