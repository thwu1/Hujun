function [ u,v,p ] = implicit_dgs( u,v,p,f,g,v1 )
% iterate v1 times 
n = size(p,1);

for k = 1:v1
    
    [ u,v ] = implicit_gs(u,v,p,f,g);
    % update iterior block
    
    for i = 2:n-1
        for j = 2:n-1
    r = u(i,j+1) - u(i,j) + v(i+1,j) - v(i,j);
    u(i,j+1) = u(i,j+1) - r/4;    
    u(i,j) = u(i,j) + r/4;
    v(i+1,j) = v(i+1,j) - r/4;
    v(i,j) = v(i,j) + r/4;
    r = n*r;
    p(i,j) = p(i,j) - r;
    p(i+1,j) = p(i+1,j) + r/4;
    p(i-1,j) = p(i-1,j) + r/4;
    p(i,j+1) = p(i,j+1) + r/4;
    p(i,j-1) = p(i,j-1) + r/4;
        end
    end
    % update lower boundary
    for i = 2:n-1
        r = u(1,i+1) - u(1,i) + v(2,i) - v(1,i);
    u(1,i+1) = u(1,i+1) - r/3;
    u(1,i) = u(1,i) + r/3;
    v(2,i) = v(2,i) + r/3;
        r = n*r;
    p(1,i) = p(1,i) - r;
    p(2,i) = p(2,i) + r/3;
    p(1,i-1) = p(1,i-1) + r/3;
    p(1,i+1) = p(1,i+1) + r/3;
    end
    % update upper boundary
    for i = 2:n-1
    r = u(n,i+1) - u(n,i) + v(n+1,i) - v(n,i);
    u(n,i) = u(n,i) + r/3;
    u(n,i+1) = v(n,i+1) - r/3;
    v(n-1,i) = v(n-1,i) + r/3;
    r = n*r;
    p(n,i) = p(n,i) - r;
    p(n-1,i) = p(n-1,i) + r/3;
    p(n,i+1) = p(n,i+1) + r/3;
    p(n,i-1) = p(n,i-1) + r/3;
    end
    % update left boundary
    for i = 2:n-1
        r = u(i,2) - u(i,1) + v(i+1,1) - v(i,1);
    v(i,1) = v(i,1) + r/3;
    v(i+1,1) = v(i+1,1) - r/3;
    u(i,2) = u(i,2) - r/3;
        r = n*r;
    p(i,1) = p(i,1) - r;
    p(i,2) = p(i,2) + r/3;
    p(i-1,1) = p(i-1,1) + r/3;
    p(i+1,1) = p(i+1,1) + r/3;
    end

    % update right boundary
    for i = 2:n-1
        r = u(i,n+1) - u(i,n) + v(i+1,n) - v(i,n);
    v(i,n) = v(i,n) + r/3;
    v(i+1,n) = v(i+1,n) - r/3;
    u(i,n-1) = u(i,n-1) + r/3;
        r = n*r;
    p(i,n) = p(i,n) - r;
    p(i,n-1) = p(i,n-1) + r/3;
    p(i-1,n) = p(i-1,n) + r/3;
    p(i+1,n) = p(i+1,n) + r/3;
    end

    % left lower
    r = u(1,2) - u(1,1) + v(2,1) - v(1,1);
    u(1,2) = u(1,2) + r/2;
    v(2,1) = v(2,1) - r/2;
    r = n*r;
    p(1,1) = p(1,1) - r;
    p(1,2) = p(1,2) + r/2;
    p(2,1) = p(2,1) + r/2;
    
    % left upper
    r = u(n,2) - u(n,1) + v(n+1,1) - v(n,1);
    u(n,2) = u(n,2) - r/2;
    v(n,1) = v(n,1) + r/2;
    r = n*r;
    p(n,1) = p(n,1) - r;
    p(n,2) = p(n,2) + r/2;
    p(n-1,1) = p(n-1,1) + r/2;
    
    % right lower
    r = u(1,n+1) - u(1,n) + v(2,n) - v(1,n);
    u(1,n) = u(1,n) + r/2;
    v(2,n) = v(2,n) - r/2;
    r = n*r;
    p(1,n) = p(1,n) - r;
    p(2,n) = p(2,n) + r/2;
    p(1,n-1) = p(1,n-1) + r/2;
    
    % right upper
    r = u(n,n+1) - u(n,n) + v(n+1,n) - v(n,n);
    u(n,n) = u(n,n) + r/2;
    v(n,n) = v(n,n) + r/2;
    r = n*r;
    p(n,n) = p(n,n) - r;
    p(n-1,n) = p(n-1,n) + r/2;
    p(n,n-1) = p(n,n-1) + r/2;

end % end for
