function sinkhornTransport(a,b,K,U,lambda,stoppingCriterion="marginalDifference",p_norm=Inf,tolerance=.5e-2,maxIter= 5000,VERBOSE=0)
    # Compute the N Sinkhorn divergence for all the pairs
    # D = [d(a_1,b_1),...,d(a_N,b_N)]
    #--------------------------------------------------
    # Inputs:
    #--------------------------------------------------
    # a : is either
    #   - a d1 x 1 column vector in the probability simpolex (non-negative summing to 1). This is the [1-vs-N mode]
    #   - a d1 x N matrix, where each column vector is in the probability simplex.This is the [N x 1-vs-1 mode]
    # b : is a d2 x N matrix of N vectors in the probability simplex
    #
    # K is a d1 x d2 matrix, equal to exp(-lambda M), where M is the d1 x d2 matrix of pairwise distances between
    # bins described in a and bins in the b_1,...b_N histograms.
    # In the most simple case d_1=d_2 and M is simply a distance matrix (zero on the diagonal and such that m_ij < m_ik + m_kj
    #
    # U = K.*M is a d1 x d2 matrix, pre-stored to speed up the computation of the distances.
    #--------------------------------------------------
    # Optional Inputs:
    #--------------------------------------------------
    # Outputs:
    #--------------------------------------------------
    # D : vector of N dual-sinkhorn divergences, or upper bounds to the EMD.
    #
    # L : vector of N lower bounds to the original OT problem, a.k.a EMD. This is computed by using the dual variables of
    # the smoothed problem, which, when modified adequately, are feasible for the original (non-smoothed) OT dual problem
    #
    # u : d1 x N matrix of left scalings
    # v : d2 x N matrix of right scalings
    # The smoothed optimal transport between (a_i,b_i) can be recovered as
    # T_i = diag(u(:,i)) * K * diag(v(:,i));
    #
    # or, equivalently and substantially faster:
    # TO-DO: the julia syntax is not working now
    # T_i = broadcast(*,v(:,i)',(broadcast(*,u(:,i),K)))
    # Check this for details: https://docs.julialang.org/en/v0.6.1/manual/arrays/
    #
    # Relevant paper:
    # M. Cuturi,
    # Sinkhorn Distances : Lightspeed Computation of Optimal Transport,
    # Advances in Neural Information Processing Systems (NIPS) 26, 2013
    #
    # Optional Input Check:
    #--------------------------------------------------
    # First make sure that a and b are two dimensional array
    # TO-DO:
    #
    ## Checking the type of computation: 1-vs-N points or many pairs
    a_r1, a_c = size(a)
    b_r, b_c = size(b)


    if a_c == 1
        ONE_VS_N = true # We are computing [D(a,b_1), ... , D(a,b_N)]
    elseif a_c == b_c
        ONE_VS_N = false # We are computing [D(a_1,b_1), ... , D(a_N,b_N)]
    else
        error("The first parameter a is either a column vector in the probability simplex, or N column vectors in the probability simplex where N is size(b,2)")
    end
    ## Checking dimensionality:
    if b_c > b_r
        BIGN = true
    else
        BIGN = false
    end
    ## Small changes in the 1-vs-N case to go a bit faster.
    if ONE_VS_N # if computing 1-vs-N make sure all components of a are >0. Otherwise we can get rid of some   lines of K to go faster.
        I = find(a) # Changed from Matlab
        someZeroValues = false
        if length(I)<a_r1 # need to update some vectors and matrices if a does not have full support
            someZeroValues=true
            K = K[I,:]
            U = U[I,:]
            a = a[I]
        end
        ainvK= broadcast(/,K,a) # precomputation of this matrix saves a d1 x N Schur product at each iteration.
    end
    # display(ainvK[1:5,1:5])

    a_r = length(a)
    ## Fixed point counter
    compt=0;
    ## Initialization of
    u =  (1/a_r).*ones(a_r,b_c)./a_r #broadcast(/,ones(a_r,b_c),a_r)#
    v = zeros(b_r,b_c)
    if cmp(stoppingCriterion, "distanceRelativeDecrease")==0
        Dold = ones(1,b_c); # initialization of vector of distances.
    end
    while compt < maxIter
        if ONE_VS_N # 1-vs-N mode
            if BIGN
                u = 1./(ainvK*broadcast(/,b,(K'*u)))
            else
                u=1./(ainvK*broadcast(/,b,(u'*K)'))
            end
        else
            if BIGN
                u=broadcast(/,a,(K*broadcast(/,b,(u'*K)')))
            else
                u=broadcast(/,a,(K*broadcast(/,b,(K'*u))))
            end
        end
        compt=compt+1
        # check the stopping criterion every 20 fixed point iterations
        # or, if that's the case, before the final iteration to store the most
        # recent value for the matrix of right scaling factors v.
        if mod(compt,20)==1 || compt==maxIter
        # split computations to recover right and left scalings.
            if BIGN
                v=broadcast(/,b,(K'*u)) # main iteration of Sinkhorn's algorithm
            else
                v=broadcast(/,b,(u'*K)')
            end
            if ONE_VS_N # 1-vs-N mode
                u=1./(ainvK*v)
            else
                u=broadcast(/,a,(K*v))
            end
            # check stopping criterion
            if cmp(stoppingCriterion, "distanceRelativeDecrease")==0
                D=sum(broadcast(*,u,(U*v)))
                Criterion=norm(broadcast(/,D,Dold)-1,p_norm)
                if Criterion<tolerance || isnan(Criterion)
                    break
                end
                 Dold=D

            elseif cmp(stoppingCriterion,"marginalDifference")==0
                Criterion=norm(sum(abs.(broadcast(*,v,(K'*u))-b)),p_norm)
                if Criterion<tolerance || isnan(Criterion)
                    break
                end
            else
                error("Stopping Criterion not recognized")
            end# norm of all || . ||_1 differences between the marginal of the current solution with the actual marginals.
            compt=compt+1
            if VERBOSE>0
                display("Iteration : $compt, Criterion: $Criterion")
            end
            # if any(isnan(Criterion)), # stop all computation if a computation of one of the pairs goes wrong.
            #     error("NaN values have appeared during the fixed point iteration. This problem appears because of insufficient machine precision when processing computations with a regularization value of lambda that is too high. Try again with a reduced regularization parameter lambda or with a thresholded metric matrix M.");
            # end
        end
    end
    if cmp(stoppingCriterion,"marginalDifference")==0 # if we have been watching marginal differences, we need to compute the vector of distances.
        D=sum(u.*(U*v))
    end


    alpha = log.(u)
    beta = log.(v)
    [beta[t]==-Inf ? beta[t]=0 : beta[t]=beta[t] for t in eachindex(beta)]
    if ONE_VS_N
        L= (a'* alpha + sum(b.*beta))/lambda
    else
        [alpha[t]==-Inf ? alpha[t]=0 : alpha[t]=alpha[t] for t in eachindex(alpha)]
        L= (sum(broadcast(*,a,alpha)) + sum(broadcast(*,b,beta)))/lambda
    end
    # display(L)
    if ONE_VS_N && someZeroValues
        uu=u
        u=zeros(a_r1,b_c)
        u[I,:]=uu
    end
    return D, L, u, v
end
