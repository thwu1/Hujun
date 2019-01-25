function [ u,v,p ] = uni_dgs( u,v,p,f,g,v1 )
%以u,v,p为初值用dgs迭代v1次
n = size(p,1);
h = 1/n;
for k = 1:v1
% first step of DGS
s = 1;
% update u
while s
    if mod(s,10) == 0
        [ r_1,r_2 ] = cal_res(u,v,p,f,g);
%      fprintf("GS:%f\n",norm([r_1;r_2'],'fro'));
    if norm([r_1;r_2']/n,'fro') < 1e-2
%         fprintf("e:%f\n",cal_error(u,v,p));
        break;
    end
    end
    s = s + 1;
for i = 2:n-1
    for j = 2:n
            u(i,j) = (1/4)*( h^2 * f(i,j) - h*( p(i,j) - p(i,j-1) ) + u(i,j+1) + u(i,j-1) + u(i-1,j) + u(i+1,j) );
    end % end for
end % end for

i = 1;
for j = 2:n
            u(i,j) = (1/3)*( h^2 * f(i,j) - h*( p(i,j) - p(i,j-1) ) + u(i,j+1) + u(i,j-1) + u(i+1,j) );
end % end if

i = n;
for j = 2:n
            u(i,j) = (1/3)*( h^2 * f(i,j) - h*( p(i,j) - p(i,j-1) ) + u(i,j+1) + u(i,j-1) + u(i-1,j) );
end

% update v
for i = 2:n
    for j = 2:n-1
        
            v(i,j) = (1/4)*( h^2 * g(i,j) - h*( p(i,j) - p(i-1,j) ) + v(i,j+1) + v(i,j-1) + v(i-1,j) + v(i+1,j) );
    end % end for
end

j = 1;
for i = 2:n
            v(i,j) = (1/3)*( h^2 * g(i,j) - h*( p(i,j) - p(i-1,j) ) + v(i,j+1) + v(i-1,j) + v(i+1,j) );
end % end for

j = n;
for i = 2:n
            v(i,j) = (1/3)*( h^2 * g(i,j) - h*( p(i,j) - p(i-1,j) ) + v(i,j-1) + v(i-1,j) + v(i+1,j) );
end % end if

end % end while
    
    for i = 2:n-1
        for j = 2:n-1
    dif = u(i,j+1) - u(i,j) + v(i+1,j) - v(i,j);
    u(i,j+1) = u(i,j+1) - dif/4;    
    u(i,j) = u(i,j) + dif/4;
    v(i+1,j) = v(i+1,j) - dif/4;
    v(i,j) = v(i,j) + dif/4;
    dif = n*dif;
    p(i,j) = p(i,j) - dif;
    p(i+1,j) = p(i+1,j) + dif/4;
    p(i-1,j) = p(i-1,j) + dif/4;
    p(i,j+1) = p(i,j+1) + dif/4;
    p(i,j-1) = p(i,j-1) + dif/4;
        end
    end
    
    %更新下边界结点速度
    for i = 2:n-1
        dif = u(1,i+1) - u(1,i) + v(2,i) - v(1,i);
    u(1,i+1) = u(1,i+1) - dif/3;
    u(1,i) = u(1,i) + dif/3;
    v(2,i) = v(2,i) + dif/3;
        dif = n*dif;
    p(1,i) = p(1,i) - (4/3)*dif;
    p(2,i) = p(2,i) + dif/3;
    p(1,i-1) = p(1,i-1) + dif/3;
    p(1,i+1) = p(1,i+1) + dif/3;
    end
    
    %更新上边界结点速度
    for i = 2:n-1
        dif = u(n,i+1) - u(n,i) + v(n+1,i) - v(n,i);
    u(n,i) = u(n,i) + dif/3;
    u(n,i+1) = v(n,i+1) - dif/3;
    v(n-1,i) = v(n-1,i) + dif/3;
        dif = n*dif;
    p(n,i) = p(n,i) - (4/3)*dif;
    p(n-1,i) = p(n-1,i) + dif/3;
    p(n,i+1) = p(n,i+1) + dif/3;
    p(n,i-1) = p(n,i-1) + dif/3;
    end
    
    %更新左边界结点速度
    for i = 2:n-1
        dif = u(i,2) - u(i,1) + v(i+1,1) - v(i,1);
    v(i,1) = v(i,1) + dif/3;
    v(i+1,1) = v(i+1,1) - dif/3;
    u(i,2) = u(i,2) - dif/3;
        dif = n*dif;
    p(i,1) = p(i,1) - (4/3)*dif;
    p(i,2) = p(i,2) + dif/3;
    p(i-1,1) = p(i-1,1) + dif/3;
    p(i+1,1) = p(i+1,1) + dif/3;
    end
    
    %更新右边界结点
    for i = 2:n-1
        dif = u(i,n+1) - u(i,n) + v(i+1,n) - v(i,n);
    v(i,n) = v(i,n) + dif/3;
    v(i+1,n) = v(i+1,n) - dif/3;
    u(i,n-1) = u(i,n-1) + dif/3;
        dif = n*dif;
    p(i,n) = p(i,n) - (4/3)*dif;
    p(i,n-1) = p(i,n-1) + dif/3;
    p(i-1,n) = p(i-1,n) + dif/3;
    p(i+1,n) = p(i+1,n) + dif/3;
    end
    
    %更新左下顶点
    dif = u(1,2) - u(1,1) + v(2,1) - v(1,1);
    u(1,2) = u(1,2) + dif/2;
    v(2,1) = v(2,1) - dif/2;
    dif = n*dif;
    p(1,1) = p(1,1) - 2*dif;
    p(1,2) = p(1,2) + dif/2;
    p(2,1) = p(2,1) + dif/2;
    
    %更新左上顶点
    dif = u(n,2) - u(n,1) + v(n+1,1) - v(n,1);
    u(n,2) = u(n,2) - dif/2;
    v(n,1) = v(n,1) + dif/2;
    dif = n*dif;
    p(n,1) = p(n,1) - 2*dif;
    p(n,2) = p(n,2) + dif/2;
    p(n-1,1) = p(n-1,1) + dif/2;
    
    %更新右下顶点
    dif = u(1,n+1) - u(1,n) + v(2,n) - v(1,n);
    u(1,n) = u(1,n) + dif/2;
    v(2,n) = v(2,n) - dif/2;
    dif = n*dif;
    p(1,n) = p(1,n) - 2*dif;
    p(2,n) = p(2,n) + dif/2;
    p(1,n-1) = p(1,n-1) + dif/2;
    
    %更新右上顶点
    dif = u(n,n+1) - u(n,n) + v(n+1,n) - v(n,n);
    u(n,n) = u(n,n) + dif/2;
    v(n,n) = v(n,n) + dif/2;
    dif = n*dif;
    p(n,n) = p(n,n) - 2*dif;
    p(n-1,n) = p(n-1,n) + dif/2;
    p(n,n-1) = p(n,n-1) + dif/2;
    
end % end for
