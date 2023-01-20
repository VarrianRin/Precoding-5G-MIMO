%% ------------------------------- INIT ---------------------------------------
clear all;
addpath('functions_stats\');

V           = 3;
d           = 6;
Nrank       = 1;
SNR         = 15;
df          = 8;
B           = 8;
UserCount   = 40;
Scenario    = '3GPP_3D_UMa_NLOS';
% -- init conditions --
[CALC, GRAPH, LOAD_R, GRAPH_ALL, BUILD, CAPACITY, LLOYD] = whatt(1);
COMPRESSION = 5; 
            % 1 - ideal
            % 2 - eigenvectors      (d)
            % 3 - dft               (d)
            % 4 - freq copsression  (d df)
            % 5 - lloyd             (d df B)
% 1 - load saved correlations matrix and build precoder and plot
% 4 - load saved correlations matrix and build precoder
% 6 - no calculation, plot
% 7 - no calculation, plot every saved cdf with V
% 8 - generate correlation matrix and build precoder and plot


% -- directions to files (correlations and capacities) --
FlIdeal = ['loads\Ideal_rank=' num2str(Nrank) '_V' num2str(V) '.mat'];
FlCorWB = ['loads\WBCorrelation_V=' num2str(V) '.mat'];
FlCorSc = ['loads\ScCorrelation_V=' num2str(V) '.mat'];

if      COMPRESSION == 1
    FlCapac = ['loads\Capacity\Capacity_Rank=' num2str(Nrank) '_V=' num2str(V) '_IDEAL' '_SNR' num2str(SNR)];
elseif  COMPRESSION == 2
    FlCapac = ['loads\Capacity\Capacity_Rank=' num2str(Nrank) '_V=' num2str(V) '_d=' num2str(d) '_SNR' num2str(SNR)];
elseif  COMPRESSION == 3
    FlCapac = ['loads\Capacity\Capacity_Rank=' num2str(Nrank) '_V=' num2str(V) '_DFT_d=' num2str(d) '_SNR' num2str(SNR)];    
elseif  COMPRESSION == 4
    FlCapac = ['loads\Capacity\Capacity_Rank=' num2str(Nrank) '_V=' num2str(V) '_d=' num2str(d) '_df=' num2str(df) '_SNR' num2str(SNR)];
elseif  COMPRESSION == 5
    FlCapac = ['loads\Capacity\Capacity_Rank=' num2str(Nrank) '_V=' num2str(V) '_d=' num2str(d) '_df=' num2str(df) '_B=' num2str(B) '_SNR' num2str(SNR)];
else
    error('maximum COMPRESSION is 5')
end

filename = ['C:\Users\ksksks\ksks\proga\matlabsk\huawei\5G_MIMO\Stats\channels\' num2str(V) '_' ...
            Scenario 'Tx32Rx4_DLUL_' num2str(UserCount) '_etilt_7_dip' '.mat'];
load(filename);
H                 = H(:, :, 1 : 48, :, :);
[M, ~, Nsc, I, K] = size(H);
Nrb               = 4;
Nsc               = Nsc / Nrb;

R_sc              = zeros(M, M, Nsc, I, K);

if CALC
    %% ---------------------------- LOADING ----------------------------------
    % -- loading or generating correlation matrixes --
    if LOAD_R
        load(FlCorWB);
        load(FlCorSc);
    else
        R_wb = CorrMatrix(H, 1);

        for sc = 1 : Nsc
            R_sc(:, :, sc, :, :) = CorrMatrix(H(:, :, (sc-1) * Nrb + (1 : Nrb), :, :), 0);
        end

        save(FlCorWB, 'R_wb');
        save(FlCorSc, 'R_sc');
    end

    if BUILD
        %% --------------------------- BUILDING PRECODER -------------------------
        % -- building 1) dft or eigen basis
        %             2) building ideal precoder
        %             3) if rb analize we align phases
        %             4) if ideal precoding capacity is needed normalize
        %             5) generating real precoder

        Coeffs = UserProcess(R_sc, R_wb, Nrank, COMPRESSION, d, df, B, Nrb, LLOYD);
        W      = BaseProcess(R_wb, Coeffs, M, COMPRESSION, d);
    end

    % ------ W is final precoder ------

    %% -------------------------- CAPACITY CALC ------------------------------

    if CAPACITY
        
        Capacity = zeros(I, K);
        for k = 1 : K
            for i = 1 : I
                Capacity(i, k) = getCapacity(H(:, :, :, i, k), ...
                                             reshape(W(:, :, :, i, k), [M, Nrank, size(W, 3)]), ...
                                             Nrank, SNR);
            end
        end
        Capacity = reshape(Capacity, [I * K, 1]);
        save(FlCapac, 'Capacity');
    end

end

%% ---------------------------- GRAPH ----------------------------------------

if GRAPH 

    if COMPRESSION == 3

        load(FlCapac);
        title  = {['CDF of Capacity V = ' num2str(V) ', SNR = ' num2str(SNR)]};
        legend = {['DFT vectors d = ' num2str(d)]};

        best_graph(2, {Capacity}, {[]}, {2}, [], {'R-'}, ...
                   title, legend, {'Capacity'}, {'CDF'});
        grid on;
    
    elseif COMPRESSION == 1

        load(FlCapac);
        title  = {['CDF of Capacity V = ' num2str(V) ', SNR = ' num2str(SNR)]};
        legend = {['Ideal precoding']};

        best_graph(2, {Capacity}, {[]}, {2}, [], {'R-'}, ...
                   title, legend, {'Capacity'}, {'CDF'});
        grid on;
    
    elseif COMPRESSION == 2

        load(FlCapac);
        title  = {['CDF of Capacity V = ' num2str(V) ', SNR = ' num2str(SNR)]};
        legend = {['eigenvectors d = ' num2str(d)]};

        best_graph(2, {Capacity}, {[]}, {2}, [], {'R-'}, ...
                   title, legend, {'Capacity'}, {'CDF'});
        grid on;

    elseif COMPRESSION == 4

        load(FlCapac);
        title  = {['CDF of Capacity V = ' num2str(V) ', SNR = ' num2str(SNR)]};
        legend = {['dft and freq compression d = ' num2str(d) ' df = ' num2str(df)]};

        best_graph(2, {Capacity}, {[]}, {2}, [], {'R-'}, ...
                   title, legend, {'Capacity'}, {'CDF'});
        grid on;

    elseif COMPRESSION == 5

        load(FlCapac);
        title  = {['CDF of Capacity V = ' num2str(V) ', SNR = ' num2str(SNR)]};
        legend = {['dft, freq, lloyd compression d = ' num2str(d) ...
                   ' df = ' num2str(df) ' B = ' num2str(B)]};

        best_graph(2, {Capacity}, {[]}, {2}, [], {'R-'}, ...
                   title, legend, {'Capacity'}, {'CDF'});
        grid on;

    end

end









