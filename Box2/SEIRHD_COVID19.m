%  SEIR_COVID19   Simulate the spread of COVID-19 with an SEIRD model

%  Set constants

N = 1e5;   %  Total population
b = 0.4/N;  %  Transmission rate for mild infections (1/d) (%0.3 or 0.5)
a = 0.2;    %  Rate of progression from the exposed to infected class (1/d)
p = 0.2;    % Rate of progression to hospitalization (1/d)
g = 1/14;   %  Rate of recovery (1/d)
fu = 0.01;   %  Probability of death
Tmax = 300; %  Max. time for simulation (d)
Tint = 50;  % Time of interevention
eff = 0; % Intervention efficacy
E0  = 1;    %  Number of initially infected people

% derive metrics
R0 = b*N/p
r  = (-(a+p)+sqrt((p-a)^2+4*a*p*R0))/2
T2 = log(2)/r

%  Set colors
co=brewermap(2,'Set2');
Icolor = co(1,:);
Dcolor  = co(2,:);

%  Set initial conditions

F0 = [N-E0 E0 0 0 0 0];  % [S0 E0 I0 H0 R0 D0]

%  Solve

[tpre,Fpre] = ode45('SEIR_COVID19_eqns_v2',[0 Tint],F0,[],b,a,p,g,fu);

% Apply intervention
if Tint < Tmax
    Fint = Fpre(end,:);  % [S0 E0 I0 R0 D0]
    [tpost,Fpost] = ode45('SEIR_COVID19_eqns_v2',[0 Tmax-Tint],Fint,[],b*(1-eff),a,p,g,fu);
    t=[tpre;tpost+Tint];
    F=[Fpre;Fpost];
else
    t=tpre;
    F=Fpre;
end

%  Plot
% Prevalence
figure(5);
subplot(2,2,1)
hold off
plot(t,F(:,3)+F(:,4),'LineWidth',2,'Color',Icolor)
hold on
plot(t,F(:,6),'LineWidth',2,'Color',Dcolor)
xlabel('Time')
ylabel('Number')
%ylim([0 N])
legend('I','D','Location','East');
legend boxoff
box on

%Prevalence Log
subplot(2,2,2)
hold off
semilogy(t,F(:,3)+F(:,4),'LineWidth',2,'Color',Icolor)
hold on
semilogy(t,F(:,6),'LineWidth',2,'Color',Dcolor)
ylim([1 N])
xlabel('Time')
ylabel('Number')
legend('I','D','Location','East');
legend boxoff
box on

%Cumulative prevalence Log
subplot(2,2,3)
hold off
semilogy(t,F(:,3)+F(:,4)+F(:,5)+F(:,6),'LineWidth',2,'Color',Icolor)
hold on
semilogy(t,F(:,6),'LineWidth',2,'Color',Dcolor)
ylim([1 N])
xlabel('Time')
ylabel('Cumulative Number')
legend('I','D','Location','East');
legend boxoff
box on

%CFR
subplot(2,2,4)
hold off
semilogy(t,100*F(:,6)./(F(:,3)+F(:,4)+F(:,5)+F(:,6)),'LineWidth',2)
hold on
semilogy(t,100*fu*ones(size(t)),'LineWidth',2)
xlabel('Time')
ylabel('Percent')
ylim([0.1 100])
legend('Deaths:Cases','CFR','Location','East');
legend boxoff
box on

