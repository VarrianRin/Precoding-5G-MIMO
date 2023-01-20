function R = CorrMatrix(H, average_polar, ENERGETIC)

[M, ~, F, I, K] = size(H);
R = zeros(M, M, I, K);

for k = 1 : K
    for i = 1 : I
        for f = 1 : F
            R(:, :, i, k) = R(:, :, i, k) + conj(H(:, :, f, i, k)) * H(:, :, f, i, k).';
        end
    end
end

if average_polar
    R = R(1 : M/2, 1 : M/2, :, :) + R((M/2 + 1) : M, (M/2 + 1) : M, :, :);
end

if exist('ENERGETIC', 'var') && ENERGETIC == 1
    Energetic(R);
end

end

function Energy = Energetic(R)

[M, ~, I, K] = size(R);
Energy       = zeros(M, I, K);
S            = zeros(M, M, I, K);
U            = zeros(M, M, I, K);

for k = 1 : K
    for i = 1 : I
        [U(:, :, i, k), S(:, :, i, k), ~] = svd(R(:, :, i, k));
    end
end

for k = 1 : K
    for i = 1 : I
        summ = sum(diag(S(:, :, i, k)));

        for m = 1 : M
            Energy(m, i, k) = sum(diag(S(1 : m, 1 : m, i, k))) / summ;
        end
    end
end

% d = 1;
% OneEnergy = squeeze(Energy(d, :, :));
% figure;
% plot(sort(reshape(OneEnergy, [K * I, 1])), 0 : 1/numel(OneEnergy) : 1 - 1/numel(OneEnergy), 'LineWidth', 2);
% leg = ['Energy of first ' num2str(d) ' eigenvectors'; 'Energy of first ' '6' ' eigenvectors'; 'Energy of first ' '5' ' eigenvectors'; ...
%        'Energy of first ' '4' ' eigenvectors'; 'Energy of first ' '3' ' eigenvectors'; 'Energy of first ' '2' ' eigenvectors'];
% legend(leg);

end