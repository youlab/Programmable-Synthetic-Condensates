function plasmid_phase_separation_Fig7a()


%%
x0=[0 1e6];
tspan=[0:0.1:50];
Cdilu=[0.01 0.1 0.25 0.5 0.75 1 2.5 5 7.5 10];

for i=1:length(Cdilu)
  
    c=Cdilu(i);
    [t,y]=ode15s(@(t,x) myfun_hill_mod(t,x,c,1),tspan,x0);
    
    plot(t,y(:,2)./sum(y,2),'-','linewidth',2)
    Y(:,i)=y(:,2)./sum(y,2);
    legstr{i}=strcat('C_{dilute}=',num2str(c));
    hold on
   
end
hold off
 xlim([0 7])
    xlabel('Time','Fontsize',20)
    ylabel('Progenitor Cell Population','Fontsize',20)
    legend(legstr)
    ax = gca;
    ax.FontSize = 15; 
end

function J=myfun_hill_mod(t,x,C_dilute,n)

%% With two selection marker=>phase separation initiated

%%
% x(2): number of cells with two plasmids 
% x(1): number of cells with only one plasmid 

%% 
% assume a logistic growth
% Ntot: the bacteria carrying capacity
% mu1: specific growth rate for cells with one plasmid
% mu2: specific growth rate for cells with two plasmid
% alpha: N2 => N1+N2

Nm=1e8;
mu2=1.5;
mu1=0.1;
d1=.1;
d2=.1;

K=1;

Alpha=K/(K+C_dilute^n);

J=zeros(length(x),1);

J(1)=mu1*x(1)*(1-(x(1)+x(2))/Nm)+mu2*Alpha*x(2)*(1-(x(1)+x(2))/Nm)-d1*x(1);
J(2)=mu2*(1-Alpha)*x(2)*(1-(x(1)+x(2))/Nm)-d2*x(2);
end