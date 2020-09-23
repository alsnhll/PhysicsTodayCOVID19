%parameters
clear all

%epidemic simulation
N=1000;
Tend=200;
Ntrj=1000;
beta=0.5;
g=0.2;

%network
nettype='gammaab';
kav=10;
par2=10;
%nettype='smallwo';
%kav=10;
%par2=0.1;
%nettype='massact';
%kav=N-1;
%par2=0;
%nettype='uniform';
%kav=10;
%par2=0;
%G=networkMake(N,nettype,kav,par2);

b=beta/kav;

%for now I just made a random network
%G=rand(N,N);

%flushot=round(rand(1,10));
flushot=zeros(1,N);

%add low weight connections
%G=G+w*ones(N,N);

%run epidemic

[n_vector,total_inf_prob]=networkSIR(N,nettype,kav,par2,flushot,Tend,Ntrj,b,g);

%avg_prob_inf=mean(total_inf_prob(1:end-1))

%figure(5)

