function [centroids] = kmeans_init(X,k)

[m,n] = size(X);
valid_cols = 1:n;

dist = zeros(k,n);
compare_dist = zeros(sum(1:k),n);
prob = ones();
centroids = zeros(m,k);
invalid_cols = zeros(1,k);
count = 1;
pick_col = randi(length(valid_cols),1);                              %creates n random numbers. picks one randomly
rand_col = valid_cols(pick_col);

for i = 1:k
    invalid_cols(1,i) = rand_col;                                    %above random number becomes invalid column as it is taken
    random_centroid = X(:,rand_col);                                 %corresponding column from data matrix is drawn and made to be centroid
    centroids(:,i) = random_centroid;                                %stores the random centroid in centroids matrix.
    valid_cols = setxor(valid_cols,invalid_cols);                    %difference between the invalid col and all nums between 1:n. All remaining are valid cols
    for g = 1:i
        compare_dist(count+(g-1),:) = sum((abs(centroids(:,g)*ones(1,n) - X)).^2,1);    %finds distance^2 between centroid and columns of X
    end
    dist(i,:) = min(compare_dist(count:count+(i-1),:),[],1);
    dist(i,:) = dist(i,:)/sum(dist(i,:));                           %floors the above distances
    prob = dist(i,:);                                               %creating probability matrix
    zero_idxs = find(prob == 0);
    find_valid = setxor(zero_idxs,1:n);
    if length(find_valid) == 1
        rand_col = find_valid;
    else
        rand_col = randsample(find_valid,1,true,prob(find_valid));
    end
    count = count + i;
end
end