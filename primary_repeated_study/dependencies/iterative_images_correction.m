function [corrected_image,ph,corrected_k,bkgd,niter,...
    thresh_bkgd,thresh_phase,diagnostics] = ...
    iterative_images_correction(raw_image,thresh_bkgd,max_iterations,...
    thresh_phase,known_bkgd,refImage,truePhase)

% ITERATIVE_IMAGES_CORRECTION
% Iterative correction of MRI image artefacts caused by random
% k-space phase fluctuations.
%
% The algorithm iteratively estimates and removes phase fluctuations
% using the signal-free background as the optimisation objective.
%
% This diagnostic version records:
%   - phase-correction convergence
%   - optimisation diagnostics
%   - computational time
%   - cumulative estimated phase
%   - background fraction
%   - image-quality metrics
%   - phase-recovery metrics
%
% INPUTS
%   raw_image       : Input complex k-space data
%   thresh_bkgd     : Background threshold
%   max_iterations  : Maximum number of outer iterations
%   thresh_phase    : Phase-convergence threshold
%   known_bkgd      : Optional known background mask
%   refImage        : Reference image for image-quality metrics
%   truePhase       : Known injected phase for phase-recovery metrics
%
% OUTPUTS
%   corrected_image : Corrected image
%   ph              : Estimated cumulative phase correction
%   corrected_k     : Corrected k-space data
%   bkgd            : Final background mask
%   niter           : Number of iterations performed
%   thresh_bkgd     : Background threshold used
%   thresh_phase    : Phase threshold used
%   diagnostics     : Structure containing iteration diagnostics
%
% NOTE
%   truePhase should use the same phase convention as the estimated
%   phase correction. For the controlled simulation, this should be
%   the appropriately transformed injected phase.

%% ============================================================
% Defaults
% =============================================================

if nargin < 2
    thresh_bkgd = 0.15;
end

if nargin < 3
    max_iterations = 20;
end

if nargin < 4
    thresh_phase = 0.15;
end

if nargin < 5
    known_bkgd = [];
end

if nargin < 6
    refImage = [];
end

if nargin < 7
    truePhase = [];
end

%% ============================================================
% Required diagnostic inputs
% =============================================================

if isempty(refImage)
    error(['refImage must be provided to calculate ', ...
        'per-iteration image-quality metrics.']);
end

if isempty(truePhase)
   truePhase = linear_correction(rph);
end

%% ============================================================
% Multidimensional stacks
% =============================================================

if ndims(raw_image) > 2

    sze = size(raw_image);

    raw_image = reshape(raw_image,...
        size(raw_image,1),size(raw_image,2),[]);

    corrected_image = zeros(size(raw_image));
    corrected_k = zeros(size(raw_image));
    bkgd = false(size(raw_image));

    ph = zeros(size(raw_image,2),size(raw_image,3));
    niter = zeros(1,size(raw_image,3));

    diagnostics = cell(1,size(raw_image,3));

    for n = 1:size(raw_image,3)

        [corrected_image(:,:,n),ph(:,n),corrected_k(:,:,n),...
            bkgd(:,:,n),niter(n),~,~,diagnostics{n}] = ...
            iterative_images_correction(...
            raw_image(:,:,n),...
            thresh_bkgd,...
            max_iterations,...
            thresh_phase,...
            known_bkgd,...
            refImage(:,:,n),...
            truePhase(:,n));

    end

    corrected_image = reshape(corrected_image,sze);
    corrected_k = reshape(corrected_k,sze);
    bkgd = reshape(bkgd,sze);

    return

end

%% ============================================================
% Initial image
% =============================================================

initial_image = ifft2c(raw_image);

%% ============================================================
% Background mask
% =============================================================

if isempty(known_bkgd)

    known_bkgd = false(size(raw_image));

    bkgd = abs(initial_image) < ...
        thresh_bkgd * max(abs(initial_image(:)));

else

    bkgd = known_bkgd == 1;

end

%% ============================================================
% Initialisation
% =============================================================

ph0 = zeros(1,size(raw_image,2));

corrected_k = raw_image;

corrected_image = initial_image;

ph = zeros(1,size(raw_image,2));

%% ============================================================
% Diagnostic arrays
% =============================================================

% Iteration number
diagnostics.iteration = [];

% Convergence diagnostics
diagnostics.std_dph = [];
diagnostics.max_abs_dph = [];
diagnostics.mean_abs_dph = [];
diagnostics.rms_dph = [];

% Optimisation diagnostics
diagnostics.fval = [];
diagnostics.exitflag = [];

% Timing diagnostics
diagnostics.time = [];
diagnostics.cumulative_time = [];

% Cumulative phase diagnostics
diagnostics.phase_mean = [];
diagnostics.phase_std = [];

% Phase-recovery performance
diagnostics.PhaseRMSE = [];
diagnostics.PhaseMAE = [];
diagnostics.PhaseCorrelation = [];

% Background diagnostics
diagnostics.bkgd_fraction = [];

% Image-quality diagnostics
diagnostics.MSE = [];
diagnostics.PSNR = [];
diagnostics.SSIM = [];

% Overall convergence status
diagnostics.converged = false;
diagnostics.max_iterations_reached = false;

%% ============================================================
% Iterative correction loop
% =============================================================

for niter = 1:max_iterations
    %% --------------------------------------------------------
    % Start timer
    % ---------------------------------------------------------

    tic

    %% --------------------------------------------------------
    % Optimise phase correction
    % ---------------------------------------------------------

    [dph,corrected_k,corrected_image,fval,exitflag] = ...
        optim_image(corrected_k,bkgd,ph0);

    iteration_time = toc;

    %% --------------------------------------------------------
    % Accumulate phase correction
    % ---------------------------------------------------------

    ph = ph + dph;

    %% --------------------------------------------------------
    % Convergence diagnostics
    % ---------------------------------------------------------

    diagnostics.iteration(niter) = niter;

    diagnostics.std_dph(niter) = ...
        std(dph);

    diagnostics.max_abs_dph(niter) = ...
        max(abs(dph));

    diagnostics.mean_abs_dph(niter) = ...
        mean(abs(dph));

    diagnostics.rms_dph(niter) = ...
        sqrt(mean(dph.^2));

    %% --------------------------------------------------------
    % Optimisation diagnostics
    % ---------------------------------------------------------

    diagnostics.fval(niter) = fval;

    diagnostics.exitflag(niter) = exitflag;

    %% --------------------------------------------------------
    % Timing diagnostics
    % ---------------------------------------------------------

    diagnostics.time(niter) = iteration_time;

    diagnostics.cumulative_time(niter) = ...
        sum(diagnostics.time(1:niter));

    %% --------------------------------------------------------
    % Cumulative phase diagnostics
    % ---------------------------------------------------------

    diagnostics.phase_mean(niter) = ...
        mean(ph);

    diagnostics.phase_std(niter) = ...
        std(ph);

    %% --------------------------------------------------------
    % Phase-recovery performance
    % ---------------------------------------------------------

    phase_error = ph(:) - truePhase(:);

    diagnostics.PhaseRMSE(niter) = ...
        sqrt(mean(phase_error.^2));

    diagnostics.PhaseMAE(niter) = ...
        mean(abs(phase_error));

    if std(ph(:)) > 0 && std(truePhase(:)) > 0

        diagnostics.PhaseCorrelation(niter) = ...
            corr(ph(:),truePhase(:));

    else

        diagnostics.PhaseCorrelation(niter) = NaN;

    end

    %% --------------------------------------------------------
    % Image-quality metrics
    % ---------------------------------------------------------

    I_reference = abs(refImage);
    I_current = abs(corrected_image);

    scale = max(I_reference(:));

    if scale == 0
        error('refImage has zero maximum intensity.');
    end

    I_reference_norm = I_reference / scale;
    I_current_norm = I_current / scale;

    diagnostics.MSE(niter) = ...
        immse(I_current_norm,I_reference_norm);

    diagnostics.PSNR(niter) = ...
        psnr(I_current_norm,I_reference_norm);

    diagnostics.SSIM(niter) = ...
        ssim(I_current_norm,I_reference_norm);

    %% --------------------------------------------------------
    % Background fraction
    % ---------------------------------------------------------

    diagnostics.bkgd_fraction(niter) = ...
        mean(bkgd(:));

    %% --------------------------------------------------------
    % Display progress
    % ---------------------------------------------------------

    fprintf(['Iteration %d: std(dph)=%.6f, ', ...
        'max|dph|=%.6f, mean|dph|=%.6f, ', ...
        'RMS(dph)=%.6f, PSNR=%.4f dB, ', ...
        'SSIM=%.4f, MSE=%.6e, ', ...
        'PhaseRMSE=%.6f, PhaseMAE=%.6f, ', ...
        'PhaseCorr=%.4f, fval=%.6e, ', ...
        'exitflag=%d, time=%.2f s\n'], ...
        niter, ...
        diagnostics.std_dph(niter), ...
        diagnostics.max_abs_dph(niter), ...
        diagnostics.mean_abs_dph(niter), ...
        diagnostics.rms_dph(niter), ...
        diagnostics.PSNR(niter), ...
        diagnostics.SSIM(niter), ...
        diagnostics.MSE(niter), ...
        diagnostics.PhaseRMSE(niter), ...
        diagnostics.PhaseMAE(niter), ...
        diagnostics.PhaseCorrelation(niter), ...
        diagnostics.fval(niter), ...
        diagnostics.exitflag(niter), ...
        diagnostics.time(niter));

    %% --------------------------------------------------------
    % Original convergence criterion
    % ---------------------------------------------------------

    if std(dph) < thresh_phase

        diagnostics.converged = true;

        break

    end

    %% --------------------------------------------------------
    % Update background mask
    % ---------------------------------------------------------

    bkgd = (known_bkgd == 1) | ...
        abs(corrected_image) < ...
        thresh_bkgd * ...
        max(max(abs(corrected_image(2:end,2:end))));

    %% --------------------------------------------------------
    % Record updated background fraction
    % ---------------------------------------------------------

    diagnostics.bkgd_fraction(niter) = ...
        mean(bkgd(:));

end

%% ============================================================
% Maximum-iteration status
% =============================================================

if niter == max_iterations && ...
        ~diagnostics.converged

    diagnostics.max_iterations_reached = true;

end

%% ============================================================
% Trim unused diagnostic entries
% =============================================================

diagnostics.iteration = ...
    diagnostics.iteration(1:niter);

diagnostics.std_dph = ...
    diagnostics.std_dph(1:niter);

diagnostics.max_abs_dph = ...
    diagnostics.max_abs_dph(1:niter);

diagnostics.mean_abs_dph = ...
    diagnostics.mean_abs_dph(1:niter);

diagnostics.rms_dph = ...
    diagnostics.rms_dph(1:niter);

diagnostics.fval = ...
    diagnostics.fval(1:niter);

diagnostics.exitflag = ...
    diagnostics.exitflag(1:niter);

diagnostics.time = ...
    diagnostics.time(1:niter);

diagnostics.cumulative_time = ...
    diagnostics.cumulative_time(1:niter);

diagnostics.phase_mean = ...
    diagnostics.phase_mean(1:niter);

diagnostics.phase_std = ...
    diagnostics.phase_std(1:niter);

diagnostics.PhaseRMSE = ...
    diagnostics.PhaseRMSE(1:niter);

diagnostics.PhaseMAE = ...
    diagnostics.PhaseMAE(1:niter);

diagnostics.PhaseCorrelation = ...
    diagnostics.PhaseCorrelation(1:niter);

diagnostics.bkgd_fraction = ...
    diagnostics.bkgd_fraction(1:niter);

diagnostics.MSE = ...
    diagnostics.MSE(1:niter);

diagnostics.PSNR = ...
    diagnostics.PSNR(1:niter);

diagnostics.SSIM = ...
    diagnostics.SSIM(1:niter);

%% ============================================================
% Final diagnostic values
% =============================================================

diagnostics.final_std_dph = ...
    diagnostics.std_dph(end);

diagnostics.final_max_abs_dph = ...
    diagnostics.max_abs_dph(end);

diagnostics.final_mean_abs_dph = ...
    diagnostics.mean_abs_dph(end);

diagnostics.final_rms_dph = ...
    diagnostics.rms_dph(end);

diagnostics.final_exitflag = ...
    diagnostics.exitflag(end);

diagnostics.final_MSE = ...
    diagnostics.MSE(end);

diagnostics.final_PSNR = ...
    diagnostics.PSNR(end);

diagnostics.final_SSIM = ...
    diagnostics.SSIM(end);

diagnostics.final_bkgd_fraction = ...
    diagnostics.bkgd_fraction(end);

diagnostics.final_PhaseRMSE = ...
    diagnostics.PhaseRMSE(end);

diagnostics.final_PhaseMAE = ...
    diagnostics.PhaseMAE(end);

diagnostics.final_PhaseCorrelation = ...
    diagnostics.PhaseCorrelation(end);

end