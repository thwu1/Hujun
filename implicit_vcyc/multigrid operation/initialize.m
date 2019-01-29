function [ u,v,p ] = initialize( n )
    u = zeros(n,n+1);
    v = zeros(n+1,n);
    p = zeros(n,n);