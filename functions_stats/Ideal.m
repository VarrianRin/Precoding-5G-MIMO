function U = Ideal(R, Nrank)

[M, ~, Nrb, I, K] = size(R);
U = zeros(M, M, Nrb, I, K);

for k = 1 : K
    for i = 1 : I
        for rb = 1 : Nrb
            [U(:, :, rb, i, k), ~, ~] = svd(R(:, :, rb, i, k));
        end
    end
end

U = U(:, 1 : Nrank, :, :, :);

end