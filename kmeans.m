function [E,revised_centroids] = iterative_kmeans(X,k)

[~,n] = size(X); %capture size of X matrix
I = eye(k);
dist = zeros(k,n);
count = 1;

[centroids] = kmeans_init(X,k); %kmeans++ iniatlization

while true
    clust_num = zeros(k,1);
    
    for i = 1:n
        dist(:,i) = (vecnorm(X(:,i) - centroids)).^2; %finds distance between i-th object and centroid and stores in i-th column of kxc dist matrix.
    end
    
    [~,idx] = min(dist,[],1);   %finds minimum of each column and the idx which corresponds to its cluster label.
    
    for j = 1:k
        clust_num(j,1) = nnz(find(idx == j));  %runs through 1->k, and finds how many times they appear in idx to find number of objects in cluster
    end
    
    S = I(idx,:);
    revised_centroids = X*S*diag(1./clust_num);         %finds revised centroids
    
    if norm(revised_centroids - centroids) < 10^-10 || count == 500  %if the change is small then stop, otherwise go till count = 50 and break.
        break;
    end
    
    centroids = revised_centroids;  %if neither of above, then continue
    count  = count  + 1;
    
end

E = I(idx,:);

labels = zeros(n,1);

for i  = 1:k
    idxs = find(E(:,i) == 1);
    labels(idxs) = i;
end

end    