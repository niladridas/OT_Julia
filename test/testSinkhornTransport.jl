# Read a,b,K,U,lambda from dat file
a = readdlm("data/a.dat")
b = readdlm("data/b.dat")
K = readdlm("data/K.dat")
U = readdlm("data/U.dat")
Lam = readdlm("data/lambda.dat")
##
include("../ot/sinkhornTransport.jl")
D, L, u, v = @time sinkhornTransport(a,b,K,U,Lam,"distanceRelativeDecrease",Inf,.5e-2,50,0)
