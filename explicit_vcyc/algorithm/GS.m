function U = GS(A,B,U,P,F)
U = U + tril(A)\(F - B*P - A*U );
end