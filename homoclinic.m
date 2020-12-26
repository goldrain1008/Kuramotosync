D=NaN(length(Ainitial),1);
    for i=1:length(Ainitial)
        A=Ainitial(i);
        while A<=0.5
        alpha=pi/2-betainitial(i);
        mu=(1+A)/2;
        nu=(1-A)/2;
    F=@(t,X)[rdyn(X(1),X(2),mu,nu,alpha);psidyn(X(1),X(2),mu,nu,alpha)];
    IC=[rinitial(i),psiinitial(i)];
    [t,XX]=ode45(F,[0 1000],IC);
    endpoint=XX(end,:);
    if sqrt((endpoint(1)-1)^2)<0.005
    D(i)=A;
    break
    end
    A=A+0.001;
        end
    end
    figure(80)
    plot(betainitial,D,'g')
    xlim([0 0.25])
    ylim([0 0.5])
    xlabel('\beta')
    ylabel('A')
    hold on
    plot(hopfcurvehopf.x(4,:),hopfcurvehopf.x(3,:),'b')
    plot(saddlecurve.x(4,1:15),saddlecurve.x(3,1:15),'r')
    plot(saddlecurve.x(4,1),saddlecurve.x(3,1),'kx','LineWidth',1)  
    legend('homoclinic','Hopf','saddle-node','BTpoint','Location','SouthEast')
 


function rdot=rdyn(r,psi,mu,nu,alpha)
rdot=((1-r^2)/2)*(mu*r*cos(alpha)+nu*cos(psi-alpha));
end 

function psidot=psidyn(r,psi,mu,nu,alpha)
psidot=((1+r^2)/(2*r))*(mu*r*sin(alpha)-nu*sin(psi-alpha))-mu*sin(alpha)-nu*r*sin(psi+alpha);
end