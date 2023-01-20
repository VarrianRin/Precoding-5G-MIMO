function [ZeroIdx, NoZeroIdx] = ZeroTaps(coeffs, df, graph)

[d2, Nrank, Nrb, I, K] = size(coeffs); % d2 = d * 2


if ~exist('graph', 'var') || graph == 0 || isempty(graph)

    ZeroIdx   = zeros(Nrb - df, Nrank, I, K);
    NoZeroIdx = zeros(df, Nrank, I, K);
    for rank = 1 : Nrank
        for k = 1 : K
            for i = 1 : I

                TimeDomC = zeros(Nrb, d2);
                for t = 1 : d2
                    TimeDomC(:, t) = ifft(squeeze(coeffs(t, rank, :, i, k)));
                end

                [~, Idx] = sort(sum(abs(TimeDomC), 2), 1, 'descend');
                ZeroIdx(:, rank, i, k)   = Idx(df+1 : Nrb);
                NoZeroIdx(:, rank, i, k) = Idx(1 : df);
            end
        end
    end

else

     cdf = zeros(Nrank, Nrb, d2, I, K);


     for rank = 1 : Nrank
        for k = 1 : K
            for i = 1 : I

                TimeDomC = zeros(Nrb, d2);
                for t = 1 : d2
                    TimeDomC(:, t) = ifft(squeeze(coeffs(t, rank, :, i, k)));
                end

                [~, Idx] = sort(sum(abs(TimeDomC), 2), 1, 'descend');
                
                for r = 1 : 8

                    NoZeroIdx = Idx(1:r);
                    ZeroIdx   = Idx(r+1 : Nrb);
                    C         = zeros(Nrb, d2);
                    for t = 1 : d2
                        C(:, t)       = ifft(squeeze(coeffs(t, rank, :, i, k)));
                        C(ZeroIdx, t) = 0;
                        Recovery      = fft(C(:, t));
                        cdf(rank, r, t, i, k) = sum(abs(Recovery)) ./ sum(abs(coeffs(t, rank, :, i, k)), 3);
                    end
                end
            end
            disp(k);
        end
     end
     filename = ['stats\freq_energetic_d=' num2str(d2/2)];
     save(filename, 'cdf');

end
end


