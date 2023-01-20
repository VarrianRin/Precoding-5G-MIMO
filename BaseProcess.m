function W = BaseProcess(R_wb, coeffs, M, COMPRESSION, d)

% COMPRESSION = 1 - ideal
%               2 - eigenvectors      (d)
%               3 - dft               (d)
%               4 - freq copsression  (d df)
%               5 - lloyd             (d df B)     

[~, Nrank, Nrb, I, K] = size(coeffs);

switch COMPRESSION
    
    case 1
        W      = coeffs;

    case 2
        P      = Basis(R_wb, d);
        W      = Projection(coeffs, P, M/2, 1);

    case 3
        P      = dftBasis(R_wb, d);
        W      = Projection(coeffs, P, M/2, 1);

    case 4
        P      = dftBasis(R_wb, d);
        coeffs = Recover(coeffs, d, Nrank, Nrb, I, K);
        W      = Projection(coeffs, P, M/2, 1);

    case 5
        P      = dftBasis(R_wb, d);
        coeffs = Recover(coeffs, d, Nrank, Nrb, I, K);
        W      = Projection(coeffs, P, M/2, 1);
       
end

W = MyNormalize(W, 3);

end

function reco = Recover(coeffs, d, Nrank, Nrb, I, K)

reco = zeros(2*d, Nrank, Nrb, I, K);
for rank = 1 : Nrank
    for k = 1 : K
        for i = 1 : I
            for t = 1 : 2*d
                reco(t, rank, :, i, k) = fft(squeeze(coeffs(t, rank, :, i, k)));
            end
        end
    end
end

end
