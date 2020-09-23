%  SEIR_COVID19   Simulate the spread of COVID-19 with an SEIRD model

%  Set constants

N = 1e5;   %  Total population
b = 0.5/N;  %  Transmission rate for mild infections (1/d) (%0.3 or 0.5)
a = 0.2;    %  Rate of progression from the exposed to infected class (1/d)
p = 0.2;    % Rate of progression to hospitalization or recovery(1/d)
fp = 0.1;  % Percent of cases that will progress to hospitalization
g = 1/14;   %  Rate of recovery from hospitalization(1/d)
fu = 0.01;   %  Probability of death
Tmax = 300; %  Max. time for simulation (d)
Tint = 301;  % Time of interevention
eff = 0; % Intervention efficacy
E0  = 1;    %  Number of initially infected people

% derive metrics
R0 = b*N/p
r  = (-(a+p)+sqrt((p-a)^2+4*a*p*R0))/2
T2 = log(2)/r

%  Set colors
co=brewermap(6,'*Spectral');

%  Set initial conditions

F0 = [N-E0 E0 0 0 0 0];  % [S0 E0 I0 H0 R0 D0]

%  Solve

[tpre,Fpre] = ode45('SEIR_COVID19_eqns_v3',[0 Tint],F0,[],b,a,p,fp,g,fu);

% Apply intervention
if Tint < Tmax
    Fint = Fpre(end,:);  % [S0 E0 I0 R0 D0]
    [tpost,Fpost] = ode45('SEIR_COVID19_eqns_v3',[0 Tmax-Tint],Fint,[],b*(1-eff),a,p,fp,g,fu);
    t=[tpre;tpost+Tint];
    F=[Fpre;Fpost];
else
    t=tpre;
    F=Fpre;
end

%  Plot
% Prevalence
figure(5)
set(gcf,'DefaultAxesColorOrder',co)

subplot(2,3,1)
plot(t,100*F(:,1:6)/N,'LineWidth',2)
xlabel('Time')
ylabel('Percent')
ylim([0 100])
xlim([0 Tmax])
legend('S','E','I','H','R','D','Location','East');
legend boxoff
box on

%Prevalence Log
subplot(2,3,4)
semilogy(t,100*F(:,1:6)/N,'LineWidth',2)
ylim([100/N 100])
xlim([0 Tmax])
xlabel('Time')
ylabel('Percent')
legend('S','E','I','H','R','D','Location','East');
legend boxoff
box on

% add 70% effective social distancing

Tint = 60;  % Time of interevention
eff = 0.7; % Intervention efficacy

%  Solve

[tpre,Fpre] = ode45('SEIR_COVID19_eqns_v3',[0 Tint],F0,[],b,a,p,fp,g,fu);

% Apply intervention
if Tint < Tmax
    Fint = Fpre(end,:);  % [S0 E0 I0 R0 D0]
    [tpost,Fpost] = ode45('SEIR_COVID19_eqns_v3',[0 Tmax-Tint],Fint,[],b*(1-eff),a,p,fp,g,fu);
    t=[tpre;tpost+Tint];
    F=[Fpre;Fpost];
else
    t=tpre;
    F=Fpre;
end

%  Plot
% Prevalence

subplot(2,3,2)
plot(t,100*F(:,1:6)/N,'LineWidth',2)
xlabel('Time')
ylabel('Percent')
ylim([0 5])
xlim([0 Tmax])
legend('S','E','I','H','R','D','Location','East');
legend boxoff
box on

%Prevalence Log
subplot(2,3,5)
semilogy(t,100*F(:,1:6)/N,'LineWidth',2)
ylim([100/N 100])
xlim([0 Tmax])
xlabel('Time')
ylabel('Percent')
legend('S','E','I','H','R','D','Location','East');
legend boxoff
box on

%Isolating 80% of infectious individuals 2 day into infectious period

fq = 0.9; %fraction who will eventually be quaratined
q = 1/1; %rate of entering quaratine (/day)

Tint = 60;  % Time of interevention

%  Solve

[tpre,Fpre] = ode45('SEIR_COVID19_eqns_v3',[0 Tint],F0,[],b,a,p,fp,g,fu);

% Apply intervention
if Tint < Tmax
    Fint = [Fpre(end,:) 0 0];  % [S0 E0 I0 R0 D0 Iq0 Q0]
    [tpost,Fpost] = ode45('SEIR_COVID19_eqns_v4',[0 Tmax-Tint],Fint,[],b,a,p,fp,g,fu,fq,q);
    t=[tpre;tpost+Tint];
    Fpost2=Fpost(:,1:6);
    Fpost2(:,3)=Fpost2(:,3)+Fpost(:,7)+Fpost(:,8); %all infected individuals, quarantine or not
    F=[Fpre;Fpost2];
else
    t=tpre;
    F=Fpre;
end

%  Plot
% Prevalence

subplot(2,3,3)
plot(t,100*F(:,1:6)/N,'LineWidth',2)
xlabel('Time')
ylabel('Percent')
ylim([0 5])
xlim([0 Tmax])
legend('S','E','I','H','R','D','Location','East');
legend boxoff
box on

%Prevalence Log
subplot(2,3,6)
semilogy(t,100*F(:,1:6)/N,'LineWidth',2)
ylim([100/N 100])
xlim([0 Tmax])
xlabel('Time')
ylabel('Percent')
legend('S','E','I','H','R','D','Location','East');
legend boxoff
box on