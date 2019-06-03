# Read a,b,K,U,lambda from dat file
a = readdlm("/Users/niladridas/Documents/GITHUB/JOT/data/a.dat")
b = readdlm("/Users/niladridas/Documents/GITHUB/JOT/data/b.dat")
K = readdlm("/Users/niladridas/Documents/GITHUB/JOT/data/K.dat")
U = readdlm("/Users/niladridas/Documents/GITHUB/JOT/data/U.dat")
Lam = readdlm("/Users/niladridas/Documents/GITHUB/JOT/data/lambda.dat")
##
include("/Users/niladridas/Documents/GITHUB/JOT/ot/sinkhornTransport.jl")

D, L, u, v = @time sinkhornTransport(a,b,K,U,Lam,"distanceRelativeDecrease",Inf,.5e-2,500,1)
#
# # Profile.clear()
# # @profile sinkhornTransport(a,b,K,U,Lam,"marginalDifference",Inf,.5e-2,500,0)
#
# @code_warntype sinkhornTransport(a,b,K,U,Lam,"marginalDifference",Inf,.5e-2,500,0)
p_matOTSink = diagm(vec(u))*K*diagm(vec(v))
M = size(a)[1]
p_matOTSink = p_matOTSink.*size(a)[1]
dist_OTSink = D
x_aOTSink = zeros(M)
for i= 1:M
    x_aOTSink[i] = sum(p_matOTSink[i,j]*samples[j] for j=1:M)
end
