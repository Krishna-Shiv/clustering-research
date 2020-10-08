function [accuracy] = cluster_accuracy(E,ground_truth)

[n,m] = size(E);

cluster_membership = zeros(n,1);

for i = 1:m
    nnz_clust = find(E(:,i));
    cluster_membership(nnz_clust) = i;
end

k_perms = perms(1:m);
find_max_sum = zeros(1,factorial(m));

for i = 1:factorial(m)
    for j = 1:n
        if k_perms(i,cluster_membership(j)) == ground_truth(j)
            find_max_sum(1,i) = find_max_sum(1,i) + 1;
        end
    end
end

find_max_sum = (1/n)*find_max_sum;
[~,idx_max] = max(find_max_sum);
accuracy = find_max_sum(idx_max);

end