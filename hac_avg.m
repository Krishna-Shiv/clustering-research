function [E] = hac_avg(X,k)

[~,n] = size(X);
    
I = eye(k);
    
W = squareform(pdist(X'));
    
Z = linkage(W,'average');
    
    dendrogram(Z)
    
idx = cluster(Z,'maxclust',k);
    
E = I(idx,:);
    
end