# Discretizes k independent AR(1) processes to discrete-state Markov chains and constructs a joint Markov chain and transition matrix
# using LazyGrids,StatsFuns
# (c) Aliaksandr Zaretski, 2021
function joint_MC(ns,rhos,sigmas,scales)
    # Inputs
    # - ns     : k-by-1 vector with desired number of states for each process if ns[i] > 2, the method of Tauchen (1986) is used
    # - rhos   : k-by-1 vector of autotocorrelation coefficients
    # - sigmas : k-by-1 vector of standard deviations
    # - scales : k-by-1 vector of standard deviation multipliers defining the extreme states as in Tauchen (1986) relevant for i where ns[i] > 2
    # Outputs
    # - xs     : k-by-1 vector with ns[i]-by-1 vectors of discretized states for each process i
    # - Ps     : k-by-1 vector with ns[i]-by-ns[i] individual transition matrices for each process i
    # - X      : prod(ns)-by-k matrix containing the joint discretized grid, each column corresponds to a specific joint state
    # - P      : prod(ns)-by-prod(ns) joint transition matrix
    k=length(ns)
    xs=Vector{Vector{Float64}}(undef,k)
    Ps=Vector{Matrix{Float64}}(undef,k)
    for i=1:k
        xs[i],Ps[i]=MC_discrete(ns[i],rhos[i],sigmas[i],scales[i])
    end
    nds=ndgrid(xs...)
    X=Matrix{Float64}(undef,prod(ns),k)
    for i=1:k
        X[:,i]=nds[i][:]
    end
    P=kron(reverse(Ps)...) # reverse to be consistent with ndgrid order
    return xs,Ps,X,P
end


# Discretizes AR(1)-process to a discrete-state Markov chain
# using StatsFuns
# (c) Aliaksandr Zaretski, 2021
function MC_discrete(n,rho,sigma,scale_sd)
    if n<1
        error("n must be a natural number!")
    elseif n==1 || sigma==0.0
        x=0.0
        P=1.0
    elseif n==2
        x=[-sigma/sqrt(1.0-rho^2.0),sigma/sqrt(1.0-rho^2.0)]
        p11=(1.0+rho)/2.0
        P=[p11 1.0-p11
            1.0-p11 p11]
    else
        xmax=scale_sd*sigma/sqrt(1.0-rho^2.0) # highest possible state
        x=range(-xmax,xmax,length=n) # vector of states
        w=x[2]-x[1]  # distance between states
        P=zeros(n,n) # transition matrix
        P[:,1]=@.normcdf((x[1]+w/2.0-rho*x)/sigma)
        for j=2:(n-1)
            P[:,j]=@.normcdf((x[j]+w/2.0-rho*x)/sigma)-normcdf((x[j]-w/2.0-rho*x)/sigma)
        end
        P[:,n]=1.0-sum(P,2)
    end
    return x,P
end


# Gets state indices of discretized Markov chain based on a vector of uniform shocks
# (c) Aliaksandr Zaretski, 2021
function get_MC_ind(shocks,P,ind_0)
    # Inputs
    # - shocks : vector of U[0,1] draws
    # - P      : Markov chain transition matrix
    # - ind_0  : index of initial state
    # Outputs
    # - ind_z  : corresponding indices of MC states
    T=length(shocks)
    cumsumP=cumsum(P,dims=2)
    ind_z=Vector{Int64}(undef,T)
    ind_z[1]=findfirst(shocks[1].<=cumsumP[ind_0,:])
    for t=2:T
        ind_z[t]=findfirst(shocks[t].<=cumsumP[ind_z[t-1],:])
    end
    return ind_z
end