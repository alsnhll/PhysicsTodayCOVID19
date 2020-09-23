%Simulates a single disease on a network

%INPUT
%G: Contact network, NxN matrix of 0-1, where N is population size
%flushot: 1 or 0 depending on whether or not someone got a flu shot (can be
%           updated to be any individual level covariate affecting susceptibility)
%Tend: Maximum simulation time - will end early if infection dies out
%Ntrj: # of runs to average over

%B1: transmission rate - probability of transmission per contact per time
                        %(days)
%g1:recovery rate - probability of recover per time (/days)
%
%OUTPUT
%n_vector = [total # infected; time disease died out]
%total_inf_prob : fraction of all runs that each individual got infected

function [n_vector,total_inf_prob]=networkSIR(N,nettype,kav,par2,flushot,Tend,Ntrj,B1,g1)

G=networkMake(N,nettype,kav,par2);

%initial figure
fig=figure(5);
clf(fig)
subplot(1,3,1)
hold on
xlabel('time')
ylabel('Number infected')

tstep=0.1; %timestep for each update of variables (days)

tvec=0:tstep:((2*Tend)+tstep); % all time values

%calculate degree of each individual
deg=zeros(N,1);
deg_friend=zeros(N,1);
for i=1:N
    
    deg(i)=sum(G(i,:));
    fr_ind=find(G(i,:)); %find all nonzero entries. 
    %fr_ind=find(G(i,:)>0.1); %find all entries >0.1 (if you added
    %background contacts)
    
    if deg(i)==0
    for j=fr_ind
        deg_friend(i)=deg_friend(i)+sum(G(j,:));
    end
    deg_friend(i)=deg_friend(i)/length(fr_ind);
    end
     
end

avgdeg=mean(deg)
avgdegfriend=mean(deg_friend)

%update network
n_vector=zeros(2,Ntrj);

R0=(mean(sum(G).^2)/mean(sum(G))-1)*B1/(B1+g1) % basic reproductive ratio of disease. If less than 1, increase B1
pemerg=1-1/R0 %probability that disease causes epidemic (vs early die off)
options = optimset('Display','off');
Rinf=fsolve(@(x)(x-(1-exp(-R0*x))),0.5,options) %predicted number of 
            %individuals who will have been infected by time infinity

total_inf_prob=zeros(1,N); 

for epi=1:Ntrj %each simulated epidemic
    
    %epi
    
    t=0;
    q=1;
    IT=zeros(size(tvec)); %total # infected
    RT=zeros(size(tvec)); %total # recovered
    %ST=1-IT-RT;
    at_qss=0; %tracks whether infection has died out or not - 0 if not
    
    I=zeros(1,N);
    R=zeros(1,N);
    
    %choose initial vaccinated - if got flu shot have 50% chance of being
    %protected from the disease
    % can change or remove depending on covariates of interest
    for i=1:N
        if flushot(i)==1
            if rand(1,1)>0.5
                R(i)=1;
            end
        end
    end
    
    Rinit=R;
    RT(1)=sum(Rinit);
    
    %randomly assign node initially infected with disease, not someone who
    %is protected
    %can change to any initiation condition
    vacced=find(R==0);
    I(vacced(ceil(rand(1,1)*length(vacced))))=1;
    IT(1)=1;
    
    %repeat for each timestep
    while at_qss ==0 && t<Tend
        
        q=q+1;
        t = t + tstep;
        
        randvec=rand(1,N); % used to determine if transitions occur

        %for each individual, find the rate of a transition occuring
        total_rate=(g1*I+B1*(1-I-R).*(I*G))*tstep;
        
        %if random # is smaller than total_rate, change state (S->I, I-> R, 0->1,1->0)
        %0 if no change, 1 if change
        change=sign(sign((1-exp(-total_rate))-randvec)+1);
        %[I',R',change']
        
        %check
        if sum(change.*R)==0
        else
            [I',R',change']
            q
            error('some R individuals are supposed to change state')
        end
        
        Itemp=I;
        I=(1-change).*I+change.*(1-I-R);
        R=(1-change).*R+change.*Itemp;
        
        %check
        %no one both infected and recovered
        if isempty(find((1-I-R)<0))==0
            [I',R']
            q
            error('some individuals are both I and R')
        end
        
        
        %total infected as a function of time
        IT(q)=sum(I);
        RT(q)=sum(R);
        
        if IT(q)==0
            %infection has died out
            at_qss=1;
            n_vector(:,epi)=[RT(q)-RT(1);t];
            total_inf_prob=total_inf_prob+R-Rinit;
        end
        
    end
    
    %plot
    subplot(1,3,1)
    plot(tvec(1:q),IT(1:q),'-',tvec(1:q),RT(1:q),'-')
    
    %regenerate network
    G=networkMake(N,nettype,kav,par2);
    
end % of repeating the simulation Ntrj times

total_inf_prob=total_inf_prob/Ntrj;

subplot(1,3,1)
ylim([0 N])
legend boxoff

%distribution of final epidemic sizes
subplot(1,3,2)
histogram(n_vector(1,:))
xlabel('Epidemic final size')
ylabel('Number of occurances')


subplot(1,3,1)
legend('infected','recovered')
legend boxoff


subplot(1,3,3)
histogram(deg)
xlabel('Degree K');
ylabel('Frequency');
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
colormap(bone)
xlim([0 max(deg)+1])







end