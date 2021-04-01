function [xs,Ps,X,P]=joint_MC(ns,rhos,sigmas,scales)
% Discretizes k independent AR(1) processes to discrete-state Markov chains
% and constructs a joint Markov chain and transition matrix
% - ns     : k-by-1 vector with desired number of states for each process;
% if ns(i) > 2, the method of Tauchen (1986) is used
% - rhos   : k-by-1 vector of autotocorrelation coefficients
% - sigmas : k-by-1 vector of standard deviations
% - scales : k-by-1 vector of standard deviation multipliers defining the
% extreme states as in Tauchen (1986); relevant for i where ns(i) > 2
% - xs     : k-by-1 cell array with discretized states for each process
% - Ps     : k-by-1 cell array with individual transition matrices
% - X      : prod(ns)-by-k matrix containing the joint discretized grid;
% each row contains a specific joint state
% - P      : prod(ns)-by-prod(ns) joint transition matrix
%
% (c) Aliaksandr Zaretski, 2020

% Define objects
k=numel(ns);
xs=cell(k,1);
Ps=cell(k,1);

% Discretize each process and construct joint transition matrix
P=1;
for i=1:k
    [xs{i},Ps{i}]=MC_discrete(ns(i),rhos(i),sigmas(i),scales(i));
    P=kron(P,Ps{i});
end

% Construct joint Markov chain
g=cell(1,k);
[g{:}]=ndgrid(xs{:});
X=cell2mat(cellfun(@(x)x(:),g,'UniformOutput',false));

end


% Discretizes an AR(1)-process to a discrete-state Markov chain
function [x,P]=MC_discrete(n,rho,sigma,scale_sd)
if n<1
    error('Non-positive number of states')
elseif n==1 || sigma==0
    x=0;
    P=1;
elseif n==2
    x=[-sigma/sqrt(1-rho^2),sigma/sqrt(1-rho^2)]';
    p11=(1+rho)/2;
    P=[p11 1-p11
        1-p11 p11];
else
    xmax=scale_sd*sigma/sqrt(1-rho^2);	% highest possible state
    x=linspace(-xmax,xmax,n)';          % vector of states
    w=x(2)-x(1);                        % distance between states
    P=zeros(n);                         % transition matrix
    P(:,1)=normcdf((x(1)+w/2-rho*x)/sigma);
    for j=2:(n-1)
        P(:,j)=normcdf((x(j)+w/2-rho*x)/sigma)-normcdf((x(j)-w/2-rho*x)/sigma);
    end
    P(:,n)=1-sum(P,2);
end
end
