function [ u,v,p,f,g ] = initialize( n )
    f = zeros(n,n+1);
    g = zeros(n+1,n);
    u = zeros(n,n+1);
    v = zeros(n+1,n);
    p = zeros(n,n);