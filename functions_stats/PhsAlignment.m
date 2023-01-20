function Walign = PhsAlignment(W)

[nTX, nL, nSC, I, K] = size(W);

Walign = zeros(nTX, nL, nSC, I, K);

for k = 1 : K
    for i = 1 : I
        for idxL = 1 : nL
            for idxSC = 1 : nSC
                
                if idxSC > 1
                    w1 = Walign(:, idxL, idxSC - 1, i, k);
                else
                    w1 = W(:, idxL, 1, i, k);
                end

                w2 = W(:, idxL, idxSC, i, k);
                pl = (w2' * w1) / abs(w2' * w1);

                Walign(:, idxL, idxSC, i, k) = (w2 * pl);
            end
        end
    end
end

end