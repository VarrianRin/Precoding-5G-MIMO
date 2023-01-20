function P = Basis(R, d)

[M, ~, I, K] = size(R);
U = zeros(M, M, I, K);
P = zeros(2*M, 2*d, I, K);

for k = 1 : K
    for i = 1 : I
        [U(:, :, i, k), ~, ~] = svd(R(:, :, i, k));
    end
end

U = U(:, 1 : d, :, :);
P(1 : M, 1 : d, :, :) = U;
P((M+1) : M*2, (d+1) : 2*d, :, :) = U;

end