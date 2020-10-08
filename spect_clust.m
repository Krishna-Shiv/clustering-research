function [E,centroids,labels] = spect_clust(X,k)

[m,n] = size(X);

G_org = X'*X;
I = eye(n);

samp = (I - ones(n,n)*(1/n));

Y = X*samp;

G = Y'*Y;

g = diag(G);

e = ones(n,1);
W = -2*G + g*e' + e*g';
W1 = -W.^(1/2);

D1 = diag(W1*ones(n,1));
L1 = (D1^(-1/2))*(D1-W1)*(D1^(-1/2));
[vecs,vals] = eig(L1);
vals = real(vals);
vecs = real(vecs);
[vals,s_ind] = sort(diag(vals),'descend');
vecs = vecs(:,s_ind(1:k));
vecs = vecs';

[E,~] = kmeans(vecs,k);

centroids = zeros(m,k);
labels = zeros(n,1);

for i = 1:k
    idxs = find(E(:,i) == 1);
    labels(idxs) = i;
    centroids(:,i) = mean(X(:,idxs),2);
end

end
