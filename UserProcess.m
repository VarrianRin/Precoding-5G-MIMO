function coeffs = UserProcess(R_sc, R_wb, Nrank, COMPRESSION, d, df, B, Nrb, LLOYD)

% COMPRESSION = 1 - ideal
%               2 - eigenvectors      (d)
%               3 - dft               (d)
%               4 - freq copsression  (d df)
%               5 - lloyd             (d df B)     

ideal_W = repelem(Ideal(R_sc, Nrank), 1, 1, Nrb);

switch COMPRESSION
    
    case 1
        coeffs              = ideal_W;

    case 2
        P                   = Basis(R_wb, d);
        coeffs              = Projection(ideal_W, P, d, 0);

    case 3
        P                   = dftBasis(R_wb, d);
        coeffs              = Projection(ideal_W, P, d, 0);

    case 4
        P                   = dftBasis(R_wb, d);
        ideal_W             = PhsAlignment(ideal_W);
        coeffs              = Projection(ideal_W, P, d, 0);
        coeffs              = FreqCompression(coeffs, df);

    case 5
        P                   = dftBasis(R_wb, d);
        ideal_W             = PhsAlignment(ideal_W);
        coeffs              = Projection(ideal_W, P, d, 0);
        [coeffs, indices]   = FreqCompression(coeffs, df);
        coeffs              = Quantization(coeffs, indices, B, df, LLOYD);
       
end

end

