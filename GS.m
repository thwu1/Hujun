function [ u,v ] = GS( u,v,p,f,g )
% first step of DGS
n = size(p,1);
h = 1/n;
% update u
while 1 
    [ r_1,r_2 ] = cal_res(u,v,p,f,g);
%      fprintf("GS:%f\n",norm([r_1;r_2'],'fro'));
    if norm([r_1;r_2'],'fro') < 1e-2
        fprintf("e:%f\n",cal_error(u,v,p));
        break;
    end
    
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