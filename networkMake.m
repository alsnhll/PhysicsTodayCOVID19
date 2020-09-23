%SYNOPSIS:
%Creates a non-weighted contact network described by the class "nettype",
%average degree "kav" and another class-dependent parameter, "par2". The
%matrix G is the adjacency matrix, with each entry being 0 or 1. The option
%to output network statistics and graphs is currently commented out but is
%at the end of the code

%INPUT: 

%N:     integar, number of nodes

%nettype: string describing network class, choose from: (details within code below)
    %{
    'randomg' : random network constructed with the Erdos-Renyi/Gilbert model
    'discExp':  a random network with an exponential degree distribution  
    'uniform': a uniform random network 
    'gammaab': random network with gamma-distributed degree with tunable variance
    'univari': random network with normally-distributed degree with tunable variance
    'powerla': a scale-free (power law) network constructed using preferential
                attachment 
    'smallwo': a small-world network with a uniform degree distribution and a tunable
                re-wiring probability
    'massact' : a fully connected network
    %}

%kav:   the average degree of individuals in the network

%par2:  another class-dependent parameter

function [G]=networkMake(N,nettype,kav,par2)

%----------------create the network--------------------------------------%


%G is the matrix of connections, if Gij=1, nodes i and j are connected,
%if Gij=0, they aren't connected. G is symmetric and Gii=0

%a random network constructed with the Erdos-Renyi/Gilbert model
%kav = average degree
%par2 not used
if nettype=='randomg'
    p=kav/(N-1);
    %choose new p so still get kav despite rounding 0->1
    pp=p+(1-p)^(N-1)/(N-1);
    degs=max(1,binornd(N-1,pp,N,1));
    G=stubconnect(degs);

% a random network with an exponential degree distribution  
%kav = average degree
%par2 not used
elseif nettype=='discExp'
    nodes=(1:1:N)'; %indeces for each node
    y=log(kav/(kav-1)); %gamma param for disc exp
    degs=floor(1/y*log(N./nodes)+1); % number of connections (degree) of each node
    %+1 because don't want anyone to have degree zer
    G=stubconnect(degs);
 
% a uniform random network  
%kav = average degree
%par2 not used
elseif nettype=='uniform'
    if mod(kav,1)==0
        degs=kav*ones(N,1); %y is
        G=stubconnect(degs);
    else
        error('kav must be a natural number for a uniform network')
    end

% random network with normally-distributed degree with tunable variance.
% not a great choice of distribution because negative degrees are chosen
% which are then rounded up to 1. Gammaab is prefered
%kav = average degree
%par2 = standard deviation, [0, inf]
elseif nettype == 'univari'
    degs=max(1,round(kav*ones(N,1)+par2*randn(N,1)));
    % ensure degree of at least one for each individual
    G=stubconnect(degs);

% random network with gamma-distributed degree with tunable variance. Due
% to difficulties discretizing the gamma distribution, the output average
% degree and variance may be slightly different from the input, so may
% require  bit of tweeking - especially 
%kav = average degree
%par2 = standard deviation, [0, inf]   
elseif nettype == 'gammaab'
    if par2==0
        error('par2 cannot be 0 in a gamma, set = 0.0001 instead')
    end
    %degs=max(1,round(gamrnd((4.5/par2)^2,kav/((4.5/par2)^2),[N 1])));
    %degs=max(1,round(gamrnd((kav/par2)^2,kav/((kav/par2)^2),[N 1])));
    %degs=max(1,round(gamrnd((kav/par2)^2+0.5,kav/((kav/par2)^2+0.5),[N 1])));
    %degs=max(1,round(gamrnd(((kav*0.92)/(1.05*par2))^2,((1.05*par2)^2)/(kav*0.92),[N 1])));
    %degs=max(1,round(gamrnd(((kav+0.5)/(par2))^2,((par2)^2)/(kav+0.5),[N 1])));
    degs=round(gamrnd(((kav-1)/(par2))^2,((par2)^2)/(kav-1),[N 1]))+1;

    
    G=stubconnect(degs);

% a scale-free (power law) network constructed using the preferential
% attachment method. This algorithm always gives a network with exponent
% v=3
elseif nettype == 'powerla'
    G=scalefree(N,round(kav/2));

% a scale-free (power law) network constructed with random connections. 
% kav = exponenet
%par 2 = minimum degree
elseif nettype == 'plavari'
    if par2<0
        error('par2 cannot be less than 0 in a power law network')
    end
    if kav<2
        if kav<1
            disp('power law doesnt converge for v<1, algorithm will continue')
        end
        disp('<k> for power law infinite for v<2, algorithm will continue')
    end
    
    %IMPORTANT: kav is actually the exponent of the network here
    %v=(2*<k>-1)/(<k>-1);
    v=kav;
    m=par2;
    expAvgK=(v-1)/(v-2)*m; %will not = <k> bc of kmax
    kmax=sqrt(N); 
    umin=(m/kmax)^(v-1);
    nodes=(1:1:N)'; %indeces for each node
    unirandnums=((nodes./N)*(1-umin))+umin;
    if max(unirandnums)>1
        error('random # out of range [umin,1], over')
    else
        if min(unirandnums)<umin
            error('random # out of range [umin,1], under')
        end
    end
    degs=floor(par2./(unirandnums.^(1/(v-1))));
    G=stubconnect(degs);
    

% a small-world network with a uniform degree distribution and a tunable
% re-wiring probability, using the Santos method
%kav = degree of all individuals
%par2 = rewiring probability, [0, 1]
elseif nettype =='smallwo'
    dist=kav/2; %dist = distance out of neighbors connected to
    if par2<0
        error('must have par2>0')
    end
    if round(dist)~=dist
        error('kav must be an even integar for small world networks')
    end
    G=small_world(N,dist,par2);

% a fully-connected network  
%kav not used
%par2 not used
elseif nettype=='massact'
    %%fill in
    degs=(N-1)*ones(N,1);
    %this is very slow but not sure why
    G=ones(N,N)-eye(N);
    display('created mass action matrix, exiting')
    return;
else
    display('Not a valid argument for nettype, must be either "uniform","discExp" or "massact"')
end

    function F=stubconnect(degs)
        
        Nsize=length(degs);
        F=sparse(Nsize,Nsize);
        degsSum=sum(degs); % total number of half-connections
        
        %make sure the number of half connections is even
        if mod(degsSum,2)==0
        else
            rand1=rand(1,1);
            degs(ceil(rand1*Nsize))=degs(ceil(rand1*Nsize))+1;
            degsSum=sum(degs);
        end
        
        %create the stubs
        
        %in edge1, each node is entered a number of times equal to its degree. Thus
        %this vector lists all half connections
        edge1=zeros(degsSum,1);
        k=1;
        for i=1:Nsize %node index
            for j=1:degs(i) %degree index,
                edge1(k)=i;
                k=k+1;
            end
        end
        k=k-1;
        
        %tests
        if k==length(edge1)
        else
            disp('problem filling edge1: not expected number of edges, k=/length(edge1)');
            sprintf('k=%d,edge1 has %d elements',k, length(edge1))
        end
        
        edge1min=min(edge1);
        if edge1min==0
            disp('there is still at least one zero entry in edge1: BAD');
            return
        end
        
        
        %-----------Pair up stubs to make edges of network-------------
        q=0; % switch, turned to one when all pairs are of the type nodeX-nodeY,
        %X=/Y, or, there is only one pair left and it is nodeX-nodeX, in which case it is ignored
        f=0; %counter
        while q==0
            
            %randomly rearrange the elements in edge1
            rp=randperm(length(edge1));
            edge1=edge1(rp,:);% a random permutation of the elements in edge1
            
            %divide the elements between edge1 and edge2, and then they will be paired
            %up, edge1(i) <-> edge2(i)
            edge2=edge1(1:length(edge1)/2);
            edge1=edge1(length(edge1)/2+1:length(edge1));
            
            %test
            edge1min=min(edge1);
            if edge1min==0
                disp('there is still at least one zero entry in edge1: BAD');
                return
            end
            edge2min=min(edge2);
            if edge2min==0
                disp('there is still at least one zero entry in edge2: BAD');
                return
            end
            
            %check that there are no nodeX-nodeX pairs
            check=abs(sign(edge1-edge2)); % will be 0 if the entries in edge1 and 2 are
            %the same, and 1 otherwise
            %check=check;
            
            %fill into G all nodeX-nodeY pairs
            for i=1:length(edge1);
                F(edge1(i),edge2(i))=F(edge1(i),edge2(i))+1*check(i);
            end
            
            %check=check;
            %rearrange all nodeX-nodeX pairs and repeat, until we only have
            %nodeX-nodeY pairs
            if min(check)==1  %no nodeX-nodeX pairs
                q=1;
            elseif length(check)==1 % if there is only one entry left in check,
                %and it is 0, then there is no way to make a nodeX-nodeY pair,
                %so just ignore this connection
                q=1;
            else % there are some nodeX-nodeX pairs
                check=abs(check-1); %will be 1 if the entries in 1 and 2 are the same, and 0 otherwise
                edge1=edge1.*check;  %'erase' the good nodeX-nodeY pairs, already entered in G
                edge2=edge2.*check;
                if edge1==edge2; else disp('edge1=/edge2, BAD');end
                %order array so the zero entries are first
                [check,IX]=sort(check);
                edge1=edge1(IX,:);
                edge2=edge2(IX,:);
                %get rid of zero entries
                edge1=edge1((length(edge1)-sum(check)+1):length(edge1));
                edge2=edge2((length(edge2)-sum(check)+1):length(edge2));
                check=check((length(check)-sum(check)+1):length(check));
                if edge1==edge2; else disp('edge1=/edge2, BAD');end
                edge1=cat(1,edge1,edge2);
            end
            
            %sprintf('sum(sum(F))=%d',full(2*sum(sum(F))))
            %sprintf('q=%d',q)
            f=f+1;
            if f>10
                disp('f>10, exiting')
                return
            end
        end
        
        F=F+F'; %reflects entries, making G symmetrical
        F=sign(F);%incase some nodes are matched twice
        %{
        if full(sum(F))==degs'
            disp('sum(F)=degs, correct number of connections for each node');
        else
            disp('sum(F)=/degs, incorrect number of connections for each node');
            sprintf('Difference in degrees is %d', full(sum(sum(F)))-sum(degs))
            %return;
        end
        %}
    end



    function F=scalefree(Nsize,min_nodes)
        
        F=sparse(Nsize,Nsize);
        edge_vec=zeros(min_nodes*Nsize*2,2);
        
        %initially generate a fully-connected seed network of size min_nodes
        Nseed=ones(min_nodes+1,min_nodes+1)-eye(min_nodes+1,min_nodes+1);
        q=1;
        for i=1:min_nodes+1
            for j=i+1:min_nodes+1
                %Nseed(i,j)=1;
                edge_vec(q,:)=[i j];
                q=q+1;
                edge_vec(q,:)=[j i];
                q=q+1;
            end
        end
                
        F(1:min_nodes+1,1:min_nodes+1)=Nseed;
        
        for ind=(min_nodes+2):Nsize
            
            %choose which nodes they'll be connected to out of total edges
            connect_to_ind=ceil(rand(min_nodes,1)*length(nonzeros(edge_vec(:,1))));
            
            %make sure connections are all unique...if not, redo
            s=1;
            while s==1
                if size(unique(edge_vec(connect_to_ind,2)))==size(edge_vec(connect_to_ind,2))
                    s=0;
                else
                    connect_to_ind=ceil(rand(min_nodes,1)*length(nonzeros(edge_vec(:,1))));
                    %ind
                end
            end
            
            F(ind,edge_vec(connect_to_ind,2)')=1;
            
            edge_vec(q:q+min_nodes-1,:)=[ind*ones(min_nodes,1) edge_vec(connect_to_ind,2)];
            edge_vec(q+min_nodes:q+2*min_nodes-1,:)=[edge_vec(connect_to_ind,2) ind*ones(min_nodes,1)];
            
            %edge_vec(q:q+min_nodes-1,:)
            q=q+2*min_nodes;

        end
        
        F=sign(F+F'); %reflects entries, making G symmetrical
        
        nextra=(2*Nsize*min_nodes-full(sum(sum(F))))/2;
        %need to make up min_nodes*(min_nodes+1) connections, since
        %connections in Nseed were already reciprocated. May take a few
        %tries
        
        s=1;
        while s==1
            
            %choose initial nodes
            start_ind=ceil(rand(nextra,1)*Nsize);
            
            %choose which nodes they'll be connected to out of total edges
            connect_to_ind=ceil(rand(nextra,1)*length(nonzeros(edge_vec(:,1))));
            
            for i=1:nextra
                
                %make sure no self-connections, if so, regenerate
                while start_ind(i)==edge_vec(connect_to_ind(i),2)
                    connect_to_ind(i)=ceil(rand(1,1)*length(nonzeros(edge_vec(:,1))));
                end
                
                F(start_ind(i),edge_vec(connect_to_ind(i),2)')=1;
            end
            F=sign(F+F');
            edge_vec(q:q+nextra-1,:)=[start_ind edge_vec(connect_to_ind,2)];
            edge_vec(q+nextra:q+2*nextra-1,:)=[edge_vec(connect_to_ind,2) start_ind];
            q=q+2*nextra;
            
            nextra=(2*Nsize*min_nodes-full(sum(sum(F))))/2;
           
            if nextra==0 
                s=0; 
            end
            
        end
                
        %full(sum(sum(F)))
        
        %maxdeg=full(max(sum(F)))
        
    end

    function F=small_world(Nsize,dist,p)
        
        F=sparse(Nsize,Nsize);
        edge_vec=zeros(Nsize*dist,2);
        e=1;
        %connect everyone to their neighbors out to distance dist
        for i=1:dist
            for j=1:Nsize-i
                F(j,j+i)=1;
                edge_vec(e,:)=[j j+i];
                e=e+1;
            end
            for j=Nsize-i+1:Nsize
                F(j,j+i-Nsize)=1;
                edge_vec(e,:)=[j j+i-Nsize];
                e=e+1;
            end
            %{
            for j=1:i
                F(j,j-i+N)=1;
                edge_vec(e,:)=[j j-i+N];
                e=e+1;
            end
            for j=(i+1):Nsize
                F(j,j-i)=1;
                edge_vec(e,:)=[j j-i];
                e=e+1;
            end
            %}
        end
        
        %remove replicated edges, because nodes are matched twice
        %edge_vec=unique(edge_vec,'rows');
        
        %connect to random individuals with probability p
        if p>1
            error('p must be < 1')
        end
        
        Nlong_edges=round(p*2*dist*Nsize); % total # of long range connections
        
        %choose Nlong_edges entries from the edge_vec, and label these as
        %indeces to rewire
        
        rewire_ind1=unique(ceil(rand(Nlong_edges,1)*length(edge_vec)));
        %edge_vec(rewire_ind1,:)
        Nlong_edges=length(rewire_ind1);
        %random permutation of indeces to determine end node of rewire
        rewire_ind2=rewire_ind1(randperm(Nlong_edges));
        %edge_vec(rewire_ind2,2)
        
        for i=1:Nlong_edges
            %remove the old connections
            F(edge_vec(rewire_ind1(i),1),edge_vec(rewire_ind1(i),2))=0;
            %randomly rewire nodes on one side of these edges to nodes on the
            %other side
            F(edge_vec(rewire_ind1(i),1),edge_vec(rewire_ind2(i),2))=1;
        end

        
        %In case any matching to self occured, remove self entires
        for i=1:Nsize
            F(i,i)=0;
        end
        
        F=F+F'; %reflects entries, making G symmetrical
        F=sign(F);%incase some nodes are matched twice
        
    end

%Test G


if trace(G)==0;
    %disp('trace of G=0, okay');
else
    error('trace of G=/0, bad');
end



%toc
%tic

G=double(G);
%phi=trace(G^3)/(sum(sum(G^2))-trace(G^2))

%plot degree distribution
%{
kvec=full(sum(G));
degsAve=mean(kvec)%full(sum(kvec))/N
std(kvec,1)

figure(1);
%subplot(1,2,1);
binplace=0:1:max(max(kvec));
[n,xout]=hist(kvec,binplace);
bar(xout,n)
xlabel('Degree K');ylabel('P(k)');
title('Degree distribution')
text(max(kvec)*0.9,max(n)*0.9,sprintf('<k>=%4.2f',full(sum(kvec))/N))
h = findobj(gca,'Type','patch');
set(h,'FaceColor','w','EdgeColor','k')
colormap(bone)
xlim([0,max(max(kvec))+1]);
%
%}

end






