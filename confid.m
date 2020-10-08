function [C] = confid(X,centroids)

[~,n] = size(X);
[~,k] = size(centroids);

C = zeros(n,k);

for i = 1:k
    C(:,i) = (vecnorm(X - centroids(:,i)).^2)'; 
end

rsum = sum(C,2);

for j = 1:n
    C(j,:) = C(j,:)/rsum(j); 
end