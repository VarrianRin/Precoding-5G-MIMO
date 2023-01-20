function [coeffs, NoZeCoeffs] = Quantization(coeffs, indices, B, df, LLOYD)

[d2, Nrank, ~, I, K] = size(coeffs);  % d2 = d * 2
NoZeCoeffs           = zeros(d2 * df, I * K, Nrank);


% building vectors [1 x (df*d2)] x (I*K) x Nrank
% normalization 
%%%%%%%%%%%%%%%%%% CCCCCCCCCCCCCCCCCCCCCCHeck %%%%%%%%%%%%%
for rank = 1 : Nrank
    for k = 1 : K
        for i = 1 : I
            NoZeCoeffs(:, (k-1)*I + i, rank) = reshape(squeeze(coeffs(:, rank, indices(:, rank, i, k), i, k)), d2 * df, 1);
            NoZeCoeffs(:, (k-1)*I + i, rank) = NoZeCoeffs(:, (k-1)*I + i, rank) / norm(NoZeCoeffs(:, (k-1)*I + i, rank), 'fro');
        end
    end
end


% lloyd = 1 -> run algorithm
% lloyf = 0 -> use saved result
if ~exist('LLOYD', 'var') || LLOYD == 0 || isempty(LLOYD)
    
    filename = ['functions_stats\stats\lloydCB_V=3_R=' num2str(Nrank) '_B=' num2str(B) '_d=' num2str(d2 / 2) '_df=' num2str(df) '.mat'];
    load(filename);

else

    CB = zeros(d2 * df, 2^B, Nrank);
    for rank = 1 : Nrank
        CB(:, :, rank) = lloyd(squeeze(NoZeCoeffs(:, :, rank)), B, 1, d2/2, df, rank, Nrank);
    end
    filename = ['functions_stats\stats\lloydCB_V=3_R=' num2str(Nrank) '_B=' num2str(B) '_d=' num2str(d2 / 2) '_df=' num2str(df) '.mat'];
    save(filename, 'CB');

end


% quantize 
for rank = 1 : Nrank
    for ik = 1 : K * I
        NoZeCoeffs(:, ik, rank) = quantize(NoZeCoeffs(:, ik, rank), CB(:, :, rank));
    end
end

% back to coeff
for rank = 1 : Nrank
    for k = 1 : K
        for i = 1 : I
            for f = 1 : df
                coeffs(:, rank, indices(f, rank, i, k), i, k) = NoZeCoeffs((d2 * (f-1) + 1) : (d2 * f), (k-1)*I + i, rank);
            end
        end
    end
end

end

