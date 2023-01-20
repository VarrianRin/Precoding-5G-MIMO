function [Capacity_m, Capacity_s, SINR] = getCapacity(H, W, Nr, Signal2Noise)
% Input:
%         H  - channel [Ntx, Nrx, Nrb]
%         W  - precoder [Ntx, Nr, Nrb]
%         Nr - rank
%         Signal2Noise - SNR in dB
% Output:
%         Capacity_m - spectral efficeny [1x1]
%         Capacity_s - spectral effeceny
%         SINR - [Nrb x Nr]
% -----------------------------------------------------------------------------        
[~, ~, Nrb] = size(H);

% channel power normalization
P = 0;
for iSC = 1 : Nrb
    P = P + sum(diag(real(conj(H(:, :, iSC)) * H(:, :, iSC).')));
end
H  = H / sqrt(P);
Pn = (10 ^ (-0.1 * Signal2Noise)) / Nrb;


% get capacity like HQ
spef = zeros(Nrb, 1);
SINR = zeros(Nrb, Nr);

P2 = 0;
for iRB = 1 : Nrb
    
    H_est = permute(H(:, :, iRB), [2 1 3]);
    H_eff = H_est * W(:, :, iRB);
    H_eff_cov = H_eff' * H_eff;
    tmp   = H_eff_cov / Pn + eye(Nr);

    P2 = P2 + trace(real(H_eff_cov));
    spef(iRB, 1) = real(log2(det(tmp)));

    % SINR calculation
    if Nr == 1
        alpha = diag(H_eff_cov);
        SINR(iRB, 1) = alpha / Pn;
    else
        TMP = H_eff_cov + eye(Nr) * Pn;
        P = inv(TMP);
        P = P * H_eff_cov;
        alpha = abs(diag(P));
        SINR(iRB, :) = (alpha ./ (1-alpha))';
    end

end %iRB


Capacity_m = mean(spef);
Capacity_s = sum(spef);

end

    





