
load zumba_names
load zumba_net

s=zumba_net(:,1)';
t=zumba_net(:,2)';

gz=digraph(s,t);

[C,~,node_cols] = unique(zumba_names, 'stable');

figure(1);
plot(gz,'Layout','force','MarkerSize',6,'NodeCData',node_cols)
axis off