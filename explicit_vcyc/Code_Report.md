Code report:

main file:
problem1
problem2
problem3
problem4

algorithm:

dgs1.m(modified DGS)
dgs2.m(origin DGS)
GS.m(Gauss-Seidel)
inexact_uzawa.m(inexact Uzawa)
uzawa.m(exact Uzawa)

error_fun:

get_error(used to calculate the F norm of U-U_exact)

multigrid_fun:

get_F.m(calculate the F with boundary condition)
initial.m(set U,P to be 0)
mac.m(calculate A,B:for details see the report)
prol.m(calculate the prolongation operator)
rest.m(calculate the restriction operator for U)
rest_G.m(calculate the restriction operator for P)

