addpath('algorithm');
addpath('error function');
addpath('multigrid operation');
% n = 128;
% v1 = 100;
% v2 = 1;
% L = 1;
% a = 0.1;
% u = zeros(n,n+1);
% v = zeros(n+1,n);
% p = zeros(n,n);
% [u,v,p] = vcycle_uzawa(u,v,p,v1,v2,L,a,n);


n = 128;
v1 = 1;
v2 = 0;
L = 5;
u = zeros(n,n+1);
v = zeros(n+1,n);
p = zeros(n,n);
[u,v,p] = vcycle_dgs(u,v,p,v1,v2,L,n);
