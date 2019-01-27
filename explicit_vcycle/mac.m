function [A,B] = mac(N)

    H = sparse(N,N);
    T = sparse(N-1,N);
    G = sparse(N-1,N-1);
    H(1,1) = 3;
    H(1,2) = -1;
    H(N,N) = 3;
    H(N,N-1) = -1;
    for i = 2:N-1
        H(i,i-1) = -1;
        H(i,i) = 4;
        H(i,i+1) = -1;
        G(i-1,i) = -1;
        G(i,i-1) = -1;
    end
    for j = 1:N-1
        T(j,j) = -1;
        T(j,j+1) = 1;
    end
    A = blkdiag(kron(eye(N-1),H)+kron(G,eye(N)),kron(H,eye(N-1))+kron(eye(N),G));
    A = N^2*A;
    B = [kron(T,eye(N));kron(eye(N),T)];
    B = N*B;
end