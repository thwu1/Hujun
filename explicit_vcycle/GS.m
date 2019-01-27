function U = GS(A,B,U,P,F)
LA = tril(A);
U = U + LA\(F - B*P - A*U );
end