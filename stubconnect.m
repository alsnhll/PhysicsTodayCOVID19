function F=stubconnect(degs)
        
        Nsize=length(degs);
        F=zeros(Nsize,Nsize);
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