function plasmid_phase_separation_Fig7c_Para_Est()

Nm=1e8;
mu2=1.5;
mu1=0.1;
d1=.1;

x_guess=[1 1.1 1.2 1e8 1.5 0.1 0.1 1];

lb=ones(3,1)*0.01;
lb=[lb;x_guess(4:end)'*0.1];
ub=ones(3,1)*100;
ub=[ub;x_guess(4:end)'*10];


fun = @(x) ODE_main(x);

opts = optimoptions(@fmincon,'Algorithm','sqp','Display','iter');

NN=50;
  p = sobolset(length(x_guess),'Skip',1e3,'Leap',1e2);
  p = scramble(p,'MatousekAffineOwen');
  Ps = p(1:3:3*NN,:)';
  Ps=10.^(Ps.*repmat(log10(ub)-log10(lb),1,NN)+repmat(log10(lb),1,NN));
  clear p
for i=1:NN
    x0=Ps(:,i);
    [x(:,i),f(i)]=fmincon(fun,x0,[],[],[],[],lb,ub,[],opts);
end
end

function [J,f]=ODE_main(P0)
%% Experiment
Yexp=[0.243097643	0.132225914	0.055809913;
     0.147715666	0.095299335	0.00493028;
      0.096338886	0.016174907	0.002579267];

Ystd_Error=[0.022925312	0.041400912	0.018306723;
            0.036562634	0.022350699	0.002440965;
            0.005706953	0.00221215	0.002210274];
%%

x0=[0 1e6];
PP=P0(4:end);
 tspan=[0 60 90 120];

opt=odeset('RelTol',1e-2,'AbsTol',1e-2,'NonNegative',[1 2]);
for i=1:3
     PP0=[P0(i);PP];
     [t,y]=ode15s(@(t,x) myfun_hill(t,x,PP0),tspan,x0);
     Y(:,i)=y(2:end,2)./(y(2:end,1)+y(2:end,2));
end
J=sqrt(sum(sum(((Y-Yexp)./Ystd_Error).^2)))/9;

if nargout>1
     SSR=sum(sum((Yexp-Y).^2));
     Yexp_t=reshape(Yexp,9,1);
     SST=sum(sum((mean(Yexp_t)-Yexp).^2));
     R2=1-SSR/SST;
    IDP_count={'V:Y','WT','S:Y'};
   count=1;
  Tspan=[0:1:120];

  clear Y
  for i=1:3
      figure(i)
PP0=[P0(i);PP];
      [t,y]=ode15s(@(t,x) myfun_hill(t,x,PP0),Tspan,x0);
      Y(:,i)=y(:,2)./(y(:,1)+y(:,2));
      plot(t,Y(:,i))
      legstr{1}=strcat('Simulation,',IDP_count{i});    
      hold on
      errorbar(tspan(2:end),Yexp(:,i),Ystd_Error(:,i),'o')
      legstr{2}=strcat('Experiment,',IDP_count{i});      
      hold off
       legend(legstr)
       xlabel('Time','Fontsize',20)
    ylabel('Cell Population','Fontsize',20)
        ax = gca;
ax.FontSize = 15; 
   end   
   f=1;
   figure(4)
 

    
   for i=1:3
          PP0=[P0(i);PP];
        [t,y]=ode15s(@(t,x) myfun_hill(t,x,PP0),Tspan,x0);
        Y=y(:,2)./(y(:,1)+y(:,2));
        YY(:,i)=Y;
      plot(t,Y); 
      hold on
   end
   hold off
   legend(IDP_count)
end
end

function J=myfun_hill(t,x,P)

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
n=1;
% % Nm=1e8;
% % mu2=1.5;
% % mu1=0.1;
% % d1=.1;
% % d2=.1;
Cdilute=P(1);
Nm=P(2);
mu2=P(3);
mu1=P(4);
d=P(5);
K=P(6);

Alpha=K/(K+Cdilute);

J=zeros(length(x),1);

J(1)=mu1*x(1)*(1-(x(1)+x(2))/Nm)+mu2*Alpha*x(2)*(1-(x(1)+x(2))/Nm)-d*x(1);
J(2)=mu2*(1-Alpha)*x(2)*(1-(x(1)+x(2))/Nm)-d*x(2);
end