function [u,v,p] = uzawa(u,v,p,f,g,v1,a)
% uzawa smoother
n = size(p,1);
for k = 1:v1
[u,v] = implicit_gs(u,v,p,f,g);
p = p - (a*n)*( u(:,2:n+1) - u(:,1:n) + v(2:n+1,:) - v(1:n,:) );
end