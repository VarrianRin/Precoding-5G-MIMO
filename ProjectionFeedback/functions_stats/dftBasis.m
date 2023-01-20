function P = dftBasis(R, d)

[M, ~, I, K] = size(R);
dftVec       = zeros(M, d, I, K);
P            = zeros(2*M, 2*d, I, K);

[CB, N1, N2, O1, O2] = codebook(M);

if d == 1
    for k = 1 : K
        for i = 1 : I
            [~, ind]           = max(diag(abs(CB' * R(:, :, i, k) * CB)));
            dftVec(:, :, i, k) = CB(:, ind).';
        end
    end
else

    CBs = zeros(O1 * O2, N1 * N2);
    l   = 1;
    for q2 = 1 : O1
        for q1 = 1 : O2

            k  = 1;
            Id = zeros(1, N1 * N2);

            for n2 = q2 : O2 : N2*O2
                for n1 = q1 : O1 : N1*O1
                    Id(k) = N1 * O1 * (n2 - 1) + n1;
                    k     = k + 1;
                end
            end

            CBs(l, :) = sort(Id);
            l         = l + 1;
        end
    end

    for k = 1 : K
        for i = 1 : I

            Psum    = zeros(1, O1 * O2);
            BeamIdx = zeros(O1 * O2, N1 * N2);
    
            for Ii = 1 : O1 * O2
                v = CB(:, CBs(Ii, :));

                Power1          = real(diag(v' * R(:, :, i, k) * v));
                [Power1, idx1]  = sort(Power1, 'descend');

                Psum(Ii)        = sum(Power1(1:d));
                BeamIdx(Ii, :)  = idx1';
            end

            [~, Imax] = max(Psum(:));
            IdexS     = CBs(Imax, :);
            CB0       = CB(:, IdexS);
            dftVec(:, :, i, k) = CB0(:, BeamIdx(Imax, 1 : d));

        end
    end
end

P(1 : M, 1 : d, :, :)             = dftVec;
P((M+1) : M*2, (d+1) : 2*d, :, :) = dftVec;

end




