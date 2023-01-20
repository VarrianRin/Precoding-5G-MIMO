function W = MyNormalize(W, flag)

[M, Nrank, Nrb, I, K] = size(W);

if flag == 1 || flag == 3
% -- fro --
    for rank = 1 : Nrank
        for k = 1 : K
            for i = 1 : I
                for f = 1 : Nrb
                    nnorm = norm(W(:, rank, f, i, k), 'fro');
                    if nnorm ~= 0
                        W(:, rank, f, i, k) = W(:, rank, f, i, k) / nnorm;
                    end
                end
            end
        end
    end

end

if flag == 2 || flag == 3
% -- NEBF --
    for k = 1 : K
        for i = 1 : I
            Power = abs(W(:, :, :, i, k)) .^ 2;
            Power = sum(Power, 2);
            Power = sum(Power, 3);
            for iNor = 1 : numel(Power)
                W(iNor, :, :, i, k) = 1 / sqrt(M) / sqrt(Power(iNor)) * W(iNor, :, :, i, k) * sqrt(Nrb);
            end
        end
    end
end

end
