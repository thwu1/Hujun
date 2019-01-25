function [ f_out,g_out ] = restrict( f,g )
%RESTRICT restrict operator for the multigrid method
%

N = size( f,1 );

% N = size(P,1);

% P_out = 0.25 * ( P( 1:2:N,1:2:N ) + P( 1:2:N,2:2:N ) + P( 2:2:N,1:2:N ) + P( 2:2:N,2:2:N ) );

g_out = 0.5 * ( g( 1:2:N+1,1:2:N ) + g( 1:2:N+1,2:2:N ) );

f_out = 0.5 * ( f( 1:2:N,1:2:N+1 ) + f( 2:2:N,1:2:N+1 ) );

end

