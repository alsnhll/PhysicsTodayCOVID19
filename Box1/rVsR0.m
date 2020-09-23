
te=[3 3 5];
ti=[5 10 10];

t2=2:1:7;
r=log(2)./t2;

R0mat=zeros(3,length(t2));

R0mat(1,:)=(1+r*te(1)).*(1+r*ti(1));
R0mat(2,:)=(1+r*te(2)).*(1+r*ti(2));
R0mat(3,:)=(1+r*te(3)).*(1+r*ti(3));

R0gmatU=zeros(3,length(t2));
R0gmatU(1,:)=R0g(r,te(1),5,ti(1),1);
R0gmatU(2,:)=R0g(r,te(2),5,ti(2),1);
R0gmatU(3,:)=R0g(r,te(3),5,ti(3),1);


R0gmatL=zeros(3,length(t2));
R0gmatL(1,:)=R0g(r,te(1),1,ti(1),5);
R0gmatL(2,:)=R0g(r,te(2),1,ti(2),5);
R0gmatL(3,:)=R0g(r,te(3),1,ti(3),5);

figure(6)
subplot(1,2,1)
plot(r,R0mat,'.','MarkerSize',20)
hold on
plot(r,R0gmatL,'v','MarkerSize',2)
plot(r,R0gmatU,'^','MarkerSize',2)
hold off
xlabel('Observed exponential growth rate, r (days^{-1})')
ylabel('Inferred R_0')
xlim([0 0.4])
box off

subplot(1,2,2)
plot(t2,R0mat,'.','MarkerSize',20)
hold on
plot(t2,R0gmatL,'v','MarkerSize',2)
plot(t2,R0gmatU,'^','MarkerSize',2)
hold off
xlabel('Observed doubling time T_2 (days)')
ylabel('Infered R_0')
xlim([0 8])
box off


function x=R0g(r,te,ne,ti,ni)
   x=(r*ti).*((r*te/ne+1).^ne)./(1-(r*ti/ni+1).^(-ni)) ;
end