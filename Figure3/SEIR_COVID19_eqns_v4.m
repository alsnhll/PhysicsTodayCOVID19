function dF = SEIR_COVID19_eqns_v4(t,F,flag,b,a,p,fp,g,fu,fq,q)

[S,E,I,H,Iq,Q] = deal(F(1),F(2),F(3),F(4),F(7),F(8));

qe=1/(1/p - 1/q); % rate of exiting quaratine

dS  = -b*(I+Iq)*S;
dE  =  b*(I+Iq)*S - a*E;
dI =  (1-fq)*a*E - p*I;
dIq = fq*a*E - q*Iq; %individuals who will eventually enter quaratine
dQ = q*Iq- qe*Q; %quaratined individuals
dH = fp*p*I + fp*qe*Q - g*H;
dR  = (1-fp)*p*I + (1-fp)*qe*Q + (1-fu)*g*H;
dD  =  fu*g*H;

dF = [dS dE dI dH dR dD dIq dQ]';
return

%days in quaratine = infectious period - days to quarantine
% = 1/p - 1/q