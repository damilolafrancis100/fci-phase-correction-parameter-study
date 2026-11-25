clear
clc
% 1. Project folder and algorithm settings

iterN = 120;
bk = [];

% Number of independent random phase realisations
nRepeats = 30;

% 2. Load reference phantom
% ============================================================

S = load('PhantomRef.mat');

refImage = double(S.phantomRef(:,:,1,1,1,1));

% Reference k-space
k_true = fft2c(refImage);

%% ============================================================
% 3. Common image normalisation
% ============================================================

I_ref_raw = abs(refImage);

scale = max(I_ref_raw(:));

if scale == 0
    error('Reference image has zero maximum intensity.');
end

I_ref = I_ref_raw / scale;

%% ============================================================
% 4. Define 100 threshold configurations
% ============================================================

bkgd_values = [0.05 0.10 0.15 0.20];

phase_values = [0.025 0.05 0.10 0.15 0.20];

[B,P] = ndgrid(bkgd_values, phase_values);

configurations = [B(:), P(:)];

nConfigurations = size(configurations,1);

%% ============================================================
% 5. Define phase-error amplitudes
% ============================================================

phase_amplitudes_deg = [2 10 20 40 60 80];

nAmplitudes = numel(phase_amplitudes_deg);

%% ============================================================
% 6. Create experiment matrix
% ============================================================

% Columns:
% 1 = Background threshold
% 2 = Phase threshold
% 3 = Phase amplitude
% 4 = Repeat number

nRuns = nConfigurations * nAmplitudes * nRepeats;

experiments = zeros(nRuns,4);

runCounter = 0;

for repeat = 1:nRepeats

    for config = 1:nConfigurations

        for a = 1:nAmplitudes

            runCounter = runCounter + 1;

            experiments(runCounter,1) = ...
                configurations(config,1);

            experiments(runCounter,2) = ...
                configurations(config,2);

            experiments(runCounter,3) = ...
                phase_amplitudes_deg(a);

            experiments(runCounter,4) = repeat;

        end
    end
end

fprintf('Threshold configurations = %d\n', nConfigurations);
fprintf('Phase amplitudes         = %d\n', nAmplitudes);
fprintf('Repeats per condition    = %d\n', nRepeats);
fprintf('Total experiments        = %d\n\n', nRuns);

%% ============================================================
% 7. Preallocate run-level results
% ============================================================

Run = zeros(nRuns,1);
Repeat = zeros(nRuns,1);
RandomSeed = zeros(nRuns,1);

BackgroundThreshold = zeros(nRuns,1);
PhaseThreshold = zeros(nRuns,1);
PhaseAmplitude_deg = zeros(nRuns,1);

Runtime_s = zeros(nRuns,1);
Iterations = zeros(nRuns,1);

Converged = false(nRuns,1);
ExitFlag = zeros(nRuns,1);

FinalStd_dph = zeros(nRuns,1);
FinalMaxAbs_dph = zeros(nRuns,1);
FinalMeanAbs_dph = zeros(nRuns,1);
FinalRMS_dph = zeros(nRuns,1);

MSE_ghost = zeros(nRuns,1);
MSE_corrected = zeros(nRuns,1);

PSNR_ghost = zeros(nRuns,1);
PSNR_corrected = zeros(nRuns,1);

SSIM_ghost = zeros(nRuns,1);
SSIM_corrected = zeros(nRuns,1);

PSNR_gain = zeros(nRuns,1);
SSIM_gain = zeros(nRuns,1);

MSE_reduction_percent = zeros(nRuns,1);

PhaseRMSE = zeros(nRuns,1);
PhaseMAE = zeros(nRuns,1);
PhaseCorrelation = zeros(nRuns,1);

BackgroundFraction = zeros(nRuns,1);

ConditionID = strings(nRuns,1);

% Variable-length iteration-level results
iterationResults = cell(nRuns,1);

%% ============================================================
% 8. Sequential experiment loop
% ============================================================

for run = 1:nRuns
    %% --------------------------------------------------------
    % Current experiment
    % ---------------------------------------------------------

    th1 = experiments(run,1);
    th2 = experiments(run,2);

    phase_amplitude_deg = experiments(run,3);

    repeat = experiments(run,4);

    phAmplitude = deg2rad(phase_amplitude_deg);

    %% --------------------------------------------------------
    % Reproducible random seed
    % ---------------------------------------------------------

    seed = 100000 * repeat + run;

    rng(seed);

    RandomSeed(run) = seed;

    %% --------------------------------------------------------
    % Start from CLEAN k-space
    % ---------------------------------------------------------

    k = k_true;

    %% --------------------------------------------------------
    % Generate random phase errors
    % ---------------------------------------------------------

    rph = rand(size(k,2),1) * phAmplitude;

    % Transform injected phase to the same phase convention
    % used by the phase-correction algorithm
    truePhase = linear_correction(rph);
    %% --------------------------------------------------------
    % Apply phase error
    % ---------------------------------------------------------

    for line = 1:size(k,2)

        k(:,line) = k(:,line) .* exp(1i*rph(line));

    end

    %% --------------------------------------------------------
    % Progress information
    % ---------------------------------------------------------

    if mod(run,10) == 0 || run == 1

        fprintf( ...
            'Starting run %4d/%4d (%.1f%%) | Bkgd %.3f | Phase %.3f | Amp %3d deg\n', ...
            run, ...
            nRuns, ...
            100*run/nRuns, ...
            th1, ...
            th2, ...
            phase_amplitude_deg);

    end

    %% --------------------------------------------------------
    % Run correction algorithm
    % ---------------------------------------------------------

    [corrected_image,ph,corrected_k,bkgd,niter,...
        thresh_bkgd,thresh_phase,diagnostics] = ...
        iterative_images_correction( ...
        k,th1,iterN,th2,bk,refImage,truePhase);
    %% --------------------------------------------------------
    % Total algorithm runtime
    % ---------------------------------------------------------

    total_time = diagnostics.cumulative_time(end);

    %% --------------------------------------------------------
    % Store convergence diagnostics
    % ---------------------------------------------------------

    Converged(run) = diagnostics.converged;

    ExitFlag(run) = diagnostics.final_exitflag;

    FinalStd_dph(run) = diagnostics.final_std_dph;

    FinalMaxAbs_dph(run) = diagnostics.final_max_abs_dph;

    FinalMeanAbs_dph(run) = diagnostics.final_mean_abs_dph;

    FinalRMS_dph(run) = diagnostics.final_rms_dph;

    %% --------------------------------------------------------
    % Condition ID
    % ---------------------------------------------------------

    ConditionID(run) = sprintf( ...
        'B%.3f_P%.3f_A%d', ...
        th1, ...
        th2, ...
        phase_amplitude_deg);

    %% --------------------------------------------------------
    % Ghosted image
    % ---------------------------------------------------------

    ghostImage = abs(ifft2c(k));

    %% --------------------------------------------------------
    % Corrected image
    % ---------------------------------------------------------

    corrImage = abs(corrected_image);

    %% --------------------------------------------------------
    % Common normalisation
    % ---------------------------------------------------------

    I_ghost = ghostImage / scale;

    I_corr = corrImage / scale;

    %% --------------------------------------------------------
    % Image-quality metrics
    % ---------------------------------------------------------

    mse_ghost = immse(I_ghost,I_ref);

    mse_corr = immse(I_corr,I_ref);

    psnr_ghost = psnr(I_ghost,I_ref);

    psnr_corr = psnr(I_corr,I_ref);

    ssim_ghost = ssim(I_ghost,I_ref);

    ssim_corr = ssim(I_corr,I_ref);

    %% --------------------------------------------------------
    % Improvement metrics
    % ---------------------------------------------------------

    psnrGain = psnr_corr - psnr_ghost;

    ssimGain = ssim_corr - ssim_ghost;

    if mse_ghost > 0

        mseReduction = ...
            (mse_ghost - mse_corr) / mse_ghost * 100;

    else

        mseReduction = NaN;

    end

    %% --------------------------------------------------------
    % Phase recovery
    % ---------------------------------------------------------

    [truePhaseCorrected,~,~] = linear_correction(rph);

    estimatedPhase = ph(:);

    truePhaseCorrected = truePhaseCorrected(:);

    phaseError = estimatedPhase - truePhaseCorrected;

    phaseMAE = mean(abs(phaseError));

    phaseRMSE = sqrt(mean(phaseError.^2));

    if std(truePhaseCorrected) > 0 && ...
            std(estimatedPhase) > 0

        phaseCorr = corr( ...
            truePhaseCorrected, ...
            estimatedPhase);

    else

        phaseCorr = NaN;

    end

    %% --------------------------------------------------------
    % Background fraction
    % ---------------------------------------------------------

    backgroundFraction = mean(bkgd(:));

    %% --------------------------------------------------------
    % Store run-level results
    % ---------------------------------------------------------

    Run(run) = run;

    Repeat(run) = repeat;

    BackgroundThreshold(run) = th1;

    PhaseThreshold(run) = th2;

    PhaseAmplitude_deg(run) = phase_amplitude_deg;

    Runtime_s(run) = total_time;

    Iterations(run) = niter;

    MSE_ghost(run) = mse_ghost;

    MSE_corrected(run) = mse_corr;

    PSNR_ghost(run) = psnr_ghost;

    PSNR_corrected(run) = psnr_corr;

    SSIM_ghost(run) = ssim_ghost;

    SSIM_corrected(run) = ssim_corr;

    PSNR_gain(run) = psnrGain;

    SSIM_gain(run) = ssimGain;

    MSE_reduction_percent(run) = mseReduction;

    PhaseRMSE(run) = phaseRMSE;

    PhaseMAE(run) = phaseMAE;

    PhaseCorrelation(run) = phaseCorr;

    BackgroundFraction(run) = backgroundFraction;

    %% --------------------------------------------------------
    % Store iteration-level diagnostics
    % ---------------------------------------------------------

    nIterThisRun = numel(diagnostics.iteration);

    iterationResults{run} = table( ...
        repmat(run,nIterThisRun,1), ...
        repmat(repeat,nIterThisRun,1), ...
        repmat(seed,nIterThisRun,1), ...
        repmat(th1,nIterThisRun,1), ...
        repmat(th2,nIterThisRun,1), ...
        repmat(phase_amplitude_deg,nIterThisRun,1), ...
        diagnostics.iteration(:), ...
        diagnostics.std_dph(:), ...
        diagnostics.max_abs_dph(:), ...
        diagnostics.mean_abs_dph(:), ...
        diagnostics.rms_dph(:), ...
        diagnostics.fval(:), ...
        diagnostics.exitflag(:), ...
        diagnostics.MSE(:), ...
        diagnostics.PSNR(:), ...
        diagnostics.SSIM(:), ...
        diagnostics.bkgd_fraction(:), ...
        diagnostics.time(:), ...
        diagnostics.cumulative_time(:), ...
        'VariableNames', { ...
        'Run', ...
        'Repeat', ...
        'RandomSeed', ...
        'BackgroundThreshold', ...
        'PhaseThreshold', ...
        'PhaseAmplitude_deg', ...
        'Iteration', ...
        'Std_dph', ...
        'MaxAbs_dph', ...
        'MeanAbs_dph', ...
        'RMS_dph', ...
        'Fval', ...
        'ExitFlag', ...
        'MSE', ...
        'PSNR', ...
        'SSIM', ...
        'BackgroundFraction', ...
        'IterationTime_s', ...
        'CumulativeTime_s'});

end

%% ============================================================
% 9. Combine iteration-level results
% ============================================================

iterationResults = vertcat(iterationResults{:});

%% ============================================================
% 10. Calculate best image-quality iteration for each run
% ============================================================

uniqueRuns = unique(iterationResults.Run);

nUniqueRuns = numel(uniqueRuns);

BestPSNR = NaN(nRuns,1);
BestPSNRIteration = NaN(nRuns,1);

BestMSE = NaN(nRuns,1);
BestMSEIteration = NaN(nRuns,1);

BestSSIM = NaN(nRuns,1);
BestSSIMIteration = NaN(nRuns,1);

for r = 1:nUniqueRuns

    runID = uniqueRuns(r);

    idx = iterationResults.Run == runID;

    T = iterationResults(idx,:);

    [BestPSNR(runID),iPSNR] = max(T.PSNR);

    BestPSNRIteration(runID) = ...
        T.Iteration(iPSNR);

    [BestMSE(runID),iMSE] = min(T.MSE);

    BestMSEIteration(runID) = ...
        T.Iteration(iMSE);

    [BestSSIM(runID),iSSIM] = max(T.SSIM);

    BestSSIMIteration(runID) = ...
        T.Iteration(iSSIM);

end

%% ============================================================
% 11. Convergence versus best image-quality metrics
% ============================================================

ConvergenceIteration = Iterations;

DeltaBestPSNRIteration = ...
    ConvergenceIteration - BestPSNRIteration;

DeltaBestMSEIteration = ...
    ConvergenceIteration - BestMSEIteration;

DeltaBestSSIMIteration = ...
    ConvergenceIteration - BestSSIMIteration;

PSNR_loss_from_best = ...
    BestPSNR - PSNR_corrected;

MSE_increase_from_best_percent = ...
    (MSE_corrected - BestMSE) ./ BestMSE * 100;

SSIM_loss_from_best = ...
    BestSSIM - SSIM_corrected;

%% ============================================================
% 12. Create final run-level results table
% ============================================================

results = table( ...
    Run, ...
    Repeat, ...
    RandomSeed, ...
    ConditionID, ...
    BackgroundThreshold, ...
    PhaseThreshold, ...
    PhaseAmplitude_deg, ...
    Runtime_s, ...
    Iterations, ...
    Converged, ...
    ExitFlag, ...
    FinalStd_dph, ...
    FinalMaxAbs_dph, ...
    FinalMeanAbs_dph, ...
    FinalRMS_dph, ...
    MSE_ghost, ...
    MSE_corrected, ...
    PSNR_ghost, ...
    PSNR_corrected, ...
    SSIM_ghost, ...
    SSIM_corrected, ...
    PSNR_gain, ...
    SSIM_gain, ...
    MSE_reduction_percent, ...
    PhaseRMSE, ...
    PhaseMAE, ...
    PhaseCorrelation, ...
    BackgroundFraction, ...
    BestPSNR, ...
    BestPSNRIteration, ...
    BestMSE, ...
    BestMSEIteration, ...
    BestSSIM, ...
    BestSSIMIteration, ...
    ConvergenceIteration, ...
    DeltaBestPSNRIteration, ...
    DeltaBestMSEIteration, ...
    DeltaBestSSIMIteration, ...
    PSNR_loss_from_best, ...
    MSE_increase_from_best_percent, ...
    SSIM_loss_from_best);

%% ============================================================
% 13. Quality-control checks
% ============================================================

assert(height(results) == nRuns, ...
    'Number of result rows does not equal number of experiments.');

assert(height(iterationResults) >= nRuns, ...
    'Iteration-level results are unexpectedly empty.');

fprintf('\n============================================\n');
fprintf('QUALITY CONTROL\n');
fprintf('============================================\n');

fprintf('Expected runs                 : %d\n', nRuns);
fprintf('Actual run-level rows         : %d\n', height(results));
fprintf('Iteration-level rows          : %d\n', height(iterationResults));

fprintf('Missing PSNR                  : %d\n', ...
    sum(isnan(results.PSNR_corrected)));

fprintf('Missing SSIM                  : %d\n', ...
    sum(isnan(results.SSIM_corrected)));

fprintf('Missing phase RMSE            : %d\n', ...
    sum(isnan(results.PhaseRMSE)));

fprintf('Non-converged runs            : %d\n', ...
    sum(~results.Converged));

fprintf('Exitflag ~= 2                 : %d\n', ...
    sum(results.ExitFlag ~= 2));

fprintf('Runtime < 0                   : %d\n', ...
    sum(results.Runtime_s < 0));

fprintf('NaN phase correlation         : %d\n', ...
    sum(isnan(results.PhaseCorrelation)));

fprintf('Unique threshold configurations: %d\n', ...
    height(unique(results(:, ...
    {'BackgroundThreshold','PhaseThreshold'}))));

fprintf('Unique phase amplitudes       : %d\n', ...
    numel(unique(results.PhaseAmplitude_deg)));

fprintf('Runs where best PSNR != final : %d\n', ...
    sum(results.BestPSNRIteration ~= ...
    results.ConvergenceIteration));

fprintf('============================================\n');

%% ============================================================
% 14. Save numerical results
% ============================================================

save('ParameterSweep_20x6x30.mat', ...
    'results', ...
    'iterationResults', ...
    'experiments');

writetable(results, ...
    'ParameterSweep_20x6x30_RunResults.csv');

writetable(iterationResults, ...
    'ParameterSweep_20x6x30_IterationResults.csv');

fprintf('\nResults saved successfully.\n');

%% ============================================================
% 15. FIGURES
% ============================================================

fprintf('\n============================================\n');
fprintf('GENERATING FIGURES\n');
fprintf('============================================\n');

%% ------------------------------------------------------------
% Figure 1: PSNR gain heatmaps
% ------------------------------------------------------------

figure('Name','PSNR Gain Heatmaps','Color','w');

tiledlayout(2,5);

for a = 1:nAmplitudes

    nexttile;

    idx = results.PhaseAmplitude_deg == ...
        phase_amplitudes_deg(a);

    T = results(idx,:);

    Z = NaN(numel(bkgd_values),numel(phase_values));

    for i = 1:numel(bkgd_values)

        for j = 1:numel(phase_values)

            cellIdx = ...
                T.BackgroundThreshold == bkgd_values(i) & ...
                T.PhaseThreshold == phase_values(j);

            if any(cellIdx)

                Z(i,j) = T.PSNR_gain(cellIdx);

            end

        end

    end

    imagesc(phase_values,bkgd_values,Z);

    axis xy;

    xlabel('Phase threshold');

    ylabel('Background threshold');

    title(sprintf('%d°',phase_amplitudes_deg(a)));

    colorbar;

end

sgtitle('PSNR gain across threshold configurations');

%% ------------------------------------------------------------
% Figure 2: MSE reduction heatmaps
% ------------------------------------------------------------

figure('Name','MSE Reduction Heatmaps','Color','w');

tiledlayout(2,5);

for a = 1:nAmplitudes

    nexttile;

    idx = results.PhaseAmplitude_deg == ...
        phase_amplitudes_deg(a);

    T = results(idx,:);

    Z = NaN(numel(bkgd_values),numel(phase_values));

    for i = 1:numel(bkgd_values)

        for j = 1:numel(phase_values)

            cellIdx = ...
                T.BackgroundThreshold == bkgd_values(i) & ...
                T.PhaseThreshold == phase_values(j);

            if any(cellIdx)

                Z(i,j) = ...
                    T.MSE_reduction_percent(cellIdx);

            end

        end

    end

    imagesc(phase_values,bkgd_values,Z);

    axis xy;

    xlabel('Phase threshold');

    ylabel('Background threshold');

    title(sprintf('%d°',phase_amplitudes_deg(a)));

    colorbar;

end

sgtitle('MSE reduction across threshold configurations');

%% ------------------------------------------------------------
% Figure 3: Runtime heatmaps
% ------------------------------------------------------------

figure('Name','Runtime Heatmaps','Color','w');

tiledlayout(2,5);

for a = 1:nAmplitudes

    nexttile;

    idx = results.PhaseAmplitude_deg == ...
        phase_amplitudes_deg(a);

    T = results(idx,:);

   Z = NaN(numel(bkgd_values),numel(phase_values));

    for i = 1:numel(bkgd_values)

        for j = 1:numel(phase_values)

            cellIdx = ...
                T.BackgroundThreshold == bkgd_values(i) & ...
                T.PhaseThreshold == phase_values(j);

            if any(cellIdx)

                Z(i,j) = T.Runtime_s(cellIdx);

            end

        end

    end

    imagesc(phase_values,bkgd_values,Z);

    axis xy;

    xlabel('Phase threshold');

    ylabel('Background threshold');

    title(sprintf('%d°',phase_amplitudes_deg(a)));

    colorbar;

end

sgtitle('Runtime across threshold configurations');

%% ------------------------------------------------------------
% Figure 4: Number-of-iterations heatmaps
% ------------------------------------------------------------

figure('Name','Iteration Heatmaps','Color','w');

tiledlayout(2,5);

for a = 1:nAmplitudes

    nexttile;

    idx = results.PhaseAmplitude_deg == ...
        phase_amplitudes_deg(a);

    T = results(idx,:);

    Z = NaN(numel(bkgd_values),numel(phase_values));

    for i = 1:numel(bkgd_values)

        for j = 1:numel(phase_values)

            cellIdx = ...
                T.BackgroundThreshold == bkgd_values(i) & ...
                T.PhaseThreshold == phase_values(j);

            if any(cellIdx)

                Z(i,j) = T.Iterations(cellIdx);

            end

        end

    end

    imagesc(phase_values,bkgd_values,Z);

    axis xy;

    xlabel('Phase threshold');

    ylabel('Background threshold');

    title(sprintf('%d°',phase_amplitudes_deg(a)));

    colorbar;

end

sgtitle('Number of iterations across threshold configurations');

%% ------------------------------------------------------------
% Figure 5: Phase-recovery performance
% ------------------------------------------------------------

figure('Name','Phase Recovery Performance','Color','w');

tiledlayout(1,2);

nexttile;

gscatter( ...
    results.PhaseAmplitude_deg, ...
    results.PhaseRMSE, ...
    results.BackgroundThreshold);

xlabel('Phase amplitude (degrees)');
ylabel('Final phase RMSE (rad)');
title('Phase RMSE');
grid on;

nexttile;

gscatter( ...
    results.PhaseAmplitude_deg, ...
    results.PhaseCorrelation, ...
    results.BackgroundThreshold);

xlabel('Phase amplitude (degrees)');
ylabel('Phase correlation');
title('Phase correlation');
grid on;

sgtitle('Phase-recovery performance');

%% ------------------------------------------------------------
% Figure 6: Convergence versus image quality
% ------------------------------------------------------------

figure('Name','Convergence versus Image Quality','Color','w');

scatter( ...
    results.ConvergenceIteration, ...
    results.BestPSNRIteration, ...
    35, ...
    results.PSNR_loss_from_best, ...
    'filled');

hold on;

maxIterPlot = max(results.ConvergenceIteration);

plot( ...
    [1 maxIterPlot], ...
    [1 maxIterPlot], ...
    'k--');

xlabel('Convergence iteration');

ylabel('Best-PSNR iteration');

title('Numerical convergence versus image-quality optimum');

colorbar;

grid on;

hold off;

%% ------------------------------------------------------------
% Figure 7: Runtime versus iterations
% ------------------------------------------------------------

figure('Name','Runtime versus Iterations','Color','w');

scatter( ...
    results.Iterations, ...
    results.Runtime_s, ...
    35, ...
    results.PhaseAmplitude_deg, ...
    'filled');

xlabel('Number of iterations');

ylabel('Runtime (s)');

title('Runtime versus number of iterations');

colorbar;

grid on;

%% ------------------------------------------------------------
% Figure 8: Phase error versus image-quality improvement
% ------------------------------------------------------------

figure('Name','Phase Error versus Image Quality','Color','w');

scatter( ...
    results.PhaseRMSE, ...
    results.PSNR_gain, ...
    35, ...
    results.PhaseAmplitude_deg, ...
    'filled');

xlabel('Phase RMSE (rad)');

ylabel('PSNR gain (dB)');

title('Phase error versus image-quality improvement');

colorbar;

grid on;

%% ============================================================
% 16. Save figures
% ============================================================

saveas(1,'Figure1_PSNR_Gain_Heatmaps.png');
saveas(2,'Figure2_MSE_Reduction_Heatmaps.png');
saveas(3,'Figure3_Runtime_Heatmaps.png');
saveas(4,'Figure4_Iteration_Heatmaps.png');
saveas(5,'Figure5_Phase_Recovery.png');
saveas(6,'Figure6_Convergence_vs_ImageQuality.png');
saveas(7,'Figure7_Runtime_vs_Iterations.png');
saveas(8,'Figure8_PhaseError_vs_ImageQuality.png');

fprintf('\n============================================\n');
fprintf('COMPLETED\n');
fprintf('============================================\n');
fprintf('Experiments completed : %d\n',nRuns);
fprintf('Run-level rows         : %d\n',height(results));
fprintf('Iteration-level rows   : %d\n',height(iterationResults));
fprintf('Figures generated      : 8\n');
fprintf('============================================\n');