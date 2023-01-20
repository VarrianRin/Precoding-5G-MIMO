function [coeffs, NoZeroIdx] = FreqCompression(coeffs, df)

[d2, Nrank, ~, I, K] = size(coeffs);          % d2 = d * 2
[ZeroIdx, NoZeroIdx] = ZeroTaps(coeffs, df, 0);


for k = 1 : K
    for i = 1 : I
        for rank = 1 : Nrank
            for t = 1 : d2
                coeffs(t, rank, :, i, k) = ifft(squeeze(coeffs(t, rank, :, i, k)));
                coeffs(t, rank, ZeroIdx(:, rank, i, k), i, k) = 0;
            end
        end
    end
end

end

