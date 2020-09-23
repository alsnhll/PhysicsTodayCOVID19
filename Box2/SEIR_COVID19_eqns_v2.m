function dF = SEIR_COVID19_eqns_v2(t,F,flag,b,a,p,g,fu)

[S,E,I,H] = deal(F(1),F(2),F(3),F(4));

dS  = -b*I*S;
dE  =  b*I*S - a*E;
dI =  a*E     - p*I;
dH = p*I - g*H;
dR  =  (1-fu)*g*H;
dD  =  fu*g*H;

dF = [dS dE dI dH dR dD]';
return