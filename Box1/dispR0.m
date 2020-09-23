%using negative binomial

k=1000;
k2=1;
k3=0.5;
R0=3;
p=(1+R0/k)^(-1);
p2=(1+R0/k2)^(-1);
p3=(1+R0/k3)^(-1);

%distribution of secondary infections
w=nbinrnd(k,p,[1e6,1]);
x=nbinrnd(k2,p2,[1e6,1]);
y=nbinrnd(k3,p3,[1e6,1]);

allmax=max([max(w),max(x),max(y)]);

mxp=10;
wnew=w;
xnew=x;
ynew=y;
wnew(w>mxp)=mxp;
xnew(x>mxp)=mxp;
ynew(y>mxp)=mxp;

figure(3)
subplot(2,3,1)
histogram(wnew,0:1:(mxp+1),'Normalization','probability')
xlabel('# secondary infections')
ylabel('Frequency')
ylim([0 0.5])
title('k \rightarrow \infty')

subplot(2,3,2)
histogram(xnew,0:1:(mxp+1),'Normalization','probability')
xlabel('# secondary infections')
ylabel('Frequency')
ylim([0 0.5])
title('k=1')

subplot(2,3,3)
histogram(ynew,0:1:(mxp+1),'Normalization','probability')
xlabel('# secondary infections')
ylabel('Frequency')
ylim([0 0.5])
title('k=0.5')

%proportion of infections caused by top z of population
delz=0.01;
z=delz:delz:allmax;

prop_ppl=1-gamcdf(z,k,(1-p)/p)';
cdf_trans=delz*(1/R0)*cumsum(z.*gampdf(z,k,(1-p)/p))';
prop_trans=1-cdf_trans;
ind=find(prop_trans<0.8,1,'first');

subplot(2,3,4)
plot(100*prop_ppl,100*prop_trans)
hold on
plot(100*prop_ppl(ind),100*prop_trans(ind),'or')
plot([0 100*prop_ppl(ind)],100*prop_trans(ind)*[1 1],'--r')
plot(100*prop_ppl(ind)*[1 1],[0 100*prop_trans(ind)],'--r')
hold off
xlabel('% of infectious cases')
ylabel('% of transmission')

prop_ppl=1-gamcdf(z,k2,(1-p2)/p2)';
cdf_trans=delz*(1/R0)*cumsum(z.*gampdf(z,k2,(1-p2)/p2))';
prop_trans=1-cdf_trans;
ind=find(prop_trans<0.8,1,'first');

subplot(2,3,5)
plot(100*prop_ppl,100*prop_trans)
hold on
plot(100*prop_ppl(ind),100*prop_trans(ind),'or')
plot([0 100*prop_ppl(ind)],100*prop_trans(ind)*[1 1],'--r')
plot(100*prop_ppl(ind)*[1 1],[0 100*prop_trans(ind)],'--r')
hold off
xlabel('% of infectious cases')
ylabel('% of transmission')



prop_ppl=1-gamcdf(z,k3,(1-p3)/p3)';
cdf_trans=delz*(1/R0)*cumsum(z.*gampdf(z,k3,(1-p3)/p3))';
prop_trans=1-cdf_trans;
ind=find(prop_trans<0.8,1,'first');

subplot(2,3,6)
plot(100*prop_ppl,100*prop_trans)
hold on
plot(100*prop_ppl(ind),100*prop_trans(ind),'or')
plot([0 100*prop_ppl(ind)],100*prop_trans(ind)*[1 1],'--r')
plot(100*prop_ppl(ind)*[1 1],[0 100*prop_trans(ind)],'--r')
hold off
xlabel('% of infectious cases')
ylabel('% of transmission')

%print(gcf, '-dpdf', 'dispR0v2.pdf','-fillpage');
