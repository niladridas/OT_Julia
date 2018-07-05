function sinkhornTransport(r,c,K,U,Lam,stoppingCriterion="marginalDifference",p_norm=Inf,tolerance=.5e-2,maxIter= 5000,VERBOSE=0)
    I = find(r)
    K = K[I,:]
    U = U[I,:]
    r = r[I]
    ainvK = K./r
    compt=0;
    u = ones(size(r)[1],1)/size(r)[1]
    if stoppingCriterion == "distanceRelativeDecrease"
        Dold = 1 # initialization of vector of distances.
    end
    while compt <= maxIter
        u = 1./(ainvK*(c./(K'*u)))
        compt=compt+1;
        if mod(compt,20)==1 || compt==maxIter
            v = c./(K'*u)
            if stoppingCriterion == "distanceRelativeDecrease"
               D = sum(u.*(U*v));
               Criterion = norm(D./Dold-1,p_norm);
               if Criterion < tolerance || isnan(Criterion)
                    break;
               end
               Dold = D;
            elseif stoppingCriterion == "marginalDifference"
               D = sum(u.*(U*v));
               Criterion = norm(sum(abs(v.*(K'*u)-b)),p_norm);
               if Criterion < tolerance || isnan(Criterion)
                    break;
               end
            else
               error("Stopping Criterion not recognized");
            end
        end
    end
    return D, u, v
end
