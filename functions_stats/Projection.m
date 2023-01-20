function coeffs = Projection(ideal_W, P, d, INVERT)

[~, Nrank, Nrb, I, K] = size(ideal_W);


% -- coeffs = P' * ideal_W --
coeffs = zeros(2*d, Nrank, Nrb, I, K);
for rb = 1 : Nrb
    for k = 1 : K
        for i = 1 : I
            for rank = 1 : Nrank
                if ~INVERT
                    coeffs(:, rank, rb, i, k) = P(:, :, i, k)' * ideal_W(:, rank, rb, i, k);
                else
                    coeffs(:, rank, rb, i, k) = P(:, :, i, k) * ideal_W(:, rank, rb, i, k);
                end
            end
        end
    end
end

end