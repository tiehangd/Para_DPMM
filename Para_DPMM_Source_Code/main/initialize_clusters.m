function [clusters] = initialize_clusters(data, phi, params)

K = max(ceil(phi));  % number of clusters;

clusters = struct('logpi', [], ...
    'logpi_mn',[],'logpi_mn_l',[],'logpi_mn_r',[],...
    'logsublikelihood',[],'logsublikelihoodDelta',[],'splittable',[]);

Nk = zeros(K+1,1);
D = size(data,1);
for k=1:K
    indices = ceil(phi)==k;  % get the indexes of data points belonging to cluster k; 
    Nk(k) = nnz(indices);  % number of data points for cluster k;
    indices_l = indices & (phi-k+1<0.5); % phi seems to be one dimension for each data?
    indices_r = indices & ~indices_l;
    
    clusters(k).logpi_mn = ones(D,1)/D;
    clusters(k).logpi_mn_l = ones(D,1)/D;
    clusters(k).logpi_mn_r = ones(D,1)/D;
    clusters(k).logsublikelihood = -inf;
    clusters(k).logsublikelihoodDelta = inf;
    clusters(k).splittable = false;
end
Nk(end) = params.alpha;

[logpi] = dirrnd(Nk);
logpi = log(logpi);
for k=1:K
    clusters(k).logpi = logpi(k);
end