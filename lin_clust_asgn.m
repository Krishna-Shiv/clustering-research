function [E] = ling_clust_asgn(X,k)

beta = 0.8;
[~,n] = size(X);
I = eye(n);

samp = (I - ones(n,n)*(1/n));

Y = X*samp;

G = Y'*Y;

g = diag(G);
e = ones(n,1);

W = -2*G + g*e' + e*g';

W = W.^(1/2);

D1 = diag(W*ones(n,1));
O = (D1^(-1/2))*(W)*(D1^(-1/2));

[z,~] = svd(O,'econ');

q = (D1^(-1/2))*z(:,2);

z_k = z(:,1:k);

P = ((D1^(1/2)))*(z_k*z_k')*((D1^(1/2)));
D = diag(sqrt(1./diag(P)));

R = D*P*D;

Pc = P;
Pc(R < beta) = 0;

Pcc = Pc;
Pcc(Pc > 0) = 1;

[~,pos_q] = sort(q,'ascend');

I_rev = I(pos_q,:);

final_Pcc = I_rev*Pcc*I_rev';

m = n/k;
rho_fullstep = zeros(1,n);
rho_halfstepp = zeros(1,n);
rho_halfstepm = zeros(1,n);
count1 = 0;
count2 = 0;
count3 = 0;

for i = 1:n
    for j = 1:m
        if i - j < 1 || i + j > n
            break;
        else
            rho_fullstep(1,i) = rho_fullstep(1,i) + final_Pcc(i-j,i+j);
            count1 = count1 + 1;
        end
    end
end
for i = 1:n
    for j = 1:m
        if i - j < 1 || i + j + 1 > n
            break;
        else
            rho_halfstepp(1,i) = rho_halfstepp(1,i) + final_Pcc(i-j,i+j+1);
            count2 = count2 + 1;
        end
    end
end
for i = 1:n
    for j = 1:m
        if i - j < 1 || i + j - 1 > n
            break;
        else
            rho_halfstepm(1,i) = rho_halfstepm(1,i) + final_Pcc(i-j,i+j-1);
            count3 = count3 + 1;
        end
    end
end

a = rho_fullstep + rho_halfstepp + rho_halfstepm;

a = smooth(-a);
% a = smooth(a);

[height,peak_idxs] = findpeaks(a);

[~,idx_of_separation] = maxk(height,k-1);

actual_idxs = peak_idxs(idx_of_separation);

[f,~] = size(actual_idxs);

actual_idxs = sort(actual_idxs,'ascend');

initial = 0;
count = 1;
actual_idxs(f+1,:) = n;
[f,~] = size(actual_idxs);

trial_pos = pos_q;

labels = zeros(n,1);

for i = 1:f
    labels(initial + 1:actual_idxs(i),1) = count;
    count = count + 1;
    initial = actual_idxs(i);
    if initial == n
        break;
    end
end

labels_rev(trial_pos,1) = labels(1:end);
I_k = eye(k);
E = I_k(labels_rev,:);

end