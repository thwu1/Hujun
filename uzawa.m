function [u,v,p] = uzawa(u,v,p,f,g,v1,a)
% uzawa smoother
n = size(p,1);
for k = 1:v1
[u,v] = implicit_gs(u,v,p,f,g);
r = dif(u,v);
p = p + (a*n)*r;
end