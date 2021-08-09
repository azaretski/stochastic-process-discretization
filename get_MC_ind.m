function ind_z=get_MC_ind(shocks,P,ind_0)
% Gets state indices of discretized Markov chain based on a vector of uniform shocks
% - shocks : vector of U[0,1] draws
% - P      : Markov chain transition matrix
% - ind_0  : index of initial state
% - ind_z  : corresponding indices of MC states
%
% (c) Aliaksandr Zaretski, 2021

T=numel(shocks);
cumsumP=cumsum(P,2);
ind_z=NaN(T,1);
ind_z(1)=find(shocks(1)<=cumsumP(ind_0,:),1,'first');
for t=2:T
    ind_z(t)=find(shocks(t)<=cumsumP(ind_z(t-1),:),1,'first');
end

end