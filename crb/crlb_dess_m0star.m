  function [F, sig_M0s, sig_T1, sig_T2] = crlb_dess_m0star(M0s, T1, T2, ...
      wf, flip_spgr, flip_dess, TR_spgr, TR_dess, TE, Sig_inv)
% function [F, sig_M0s, sig_T1, sig_T2] = crlb_dess_m0star(M0s, T1, T2, ...
%     wf, flip_spgr, flip_dess, TR_spgr, TR_dess, TE, Sig_inv)
% 
% Matrix Cramer-Rao Bound calculator for unbiased estimator variance.
% This version assumes 3 unknown parameters: M0s, T1, T2.
% 
% Inputs: 
%   M0s         [1]         T2* compensated proton density, M0star
%   T1          [1]         Spin-lattice relaxation constant
%   T2          [1]         Spin-spin relaxation constant
%   wf          [1]         B0 inhomogeneity 
%   flip_spgr   [nf1]       A total nf1 of SPGR flip angle data (rad)
%   flip_dess   [nf2]       A total nf2 of DESS flip angle data (rad)
%   TR_spgr     [nf1]       SPGR repetition times (ms) 
%   TR_dess     [nf2]       DESS repetition times (ms)
%   TE          [1]         Echo time (same for both scans)
%   Sig_inv     [M M]       Inverse of noise covariance matrix
% 
% Outputs:
%   F           [3 3]       Fisher info matrix
%   sig_M0s     [1]         Time-compensated CRLB std. dev. for M0s
%   sig_T1      [1]         Time-compensated CRLB std. dev. for T1
%   sig_T2      [1]         Time-compensated CRLB std. dev. for T2
% 
% Written by: Gopal Nataraj and Jeffrey A. Fessler; Copyright 2014
% Version Two: This CRLB allow for variations in TRs between scans

% Gather constant declarations
nf1 = length(flip_spgr);
nf2 = length(flip_dess);
P = 3;                      % Number of parameters

% First, compute columns corresponding to SPGR data
grad_s = NaN(P, nf1);
for a = 1:nf1
    % Conversion to exponentials 
    E1_spgr = exp(-TR_spgr(a) ./ T1);
    E2_spgr = exp(-TR_spgr(a) ./ T2);
    
    % Append an SPGR measurement model gradient vector 
    grad_s(:,a) = ...
        [dspgr_M0s(M0s, E1_spgr, E2_spgr, flip_spgr(a), TR_spgr(a), TE);...
         dspgr_T1(M0s,  E1_spgr, E2_spgr, flip_spgr(a), TR_spgr(a), TE);...
         dspgr_T2(M0s,  E1_spgr, E2_spgr, flip_spgr(a), TR_spgr(a), TE)];
end

% Next, add columns corresponding to DESS data
grad_p = NaN(P, nf2);
grad_m = NaN(P, nf2);
for a = 1:nf2
    % Conversion to exponentials 
    E1_dess = exp(-TR_dess(a) ./ T1);
    E2_dess = exp(-TR_dess(a) ./ T2);
    
    % Append a DESS measurement model gradient vector 
    grad_p(:,a) = ...
        [dfp_M0s(M0s, E1_dess, E2_dess, flip_dess(a), TR_dess(a), TE, wf);...
         dfp_T1(M0s,  E1_dess, E2_dess, flip_dess(a), TR_dess(a), TE, wf);...
         dfp_T2(M0s,  E1_dess, E2_dess, flip_dess(a), TR_dess(a), TE, wf)];
    grad_m(:,a) = ...
        [dfm_M0s(M0s, E1_dess, E2_dess, flip_dess(a), TR_dess(a), TE, wf);...
         dfm_T1(M0s,  E1_dess, E2_dess, flip_dess(a), TR_dess(a), TE, wf);...
         dfm_T2(M0s,  E1_dess, E2_dess, flip_dess(a), TR_dess(a), TE, wf)];
end

% Fill in signal model gradient matrix, P x (M = nf1 + 2*nf2)
grad_mu = [grad_s grad_p grad_m];

% Construct Fisher information matrix 
F = grad_mu * Sig_inv * grad_mu'; 

% % Method one: time-compensated standard deviations (uses pinv())
% scan_time = sum(TR_spgr) + sum(TR_dess);
% sig_all = sqrt(scan_time) * abs(sqrt(diag(pinv(F))));
% sig_M0s = sig_all(1);
% sig_T1 = sig_all(2);
% sig_T2 = sig_all(3);

% Method two: time-compensated SDs (faster, controls cond. number)
scan_time = sum(TR_spgr) + sum(TR_dess);
max_cond_num = 1e50;
sig_M0s = sqrt(scan_time) * abs(sqrt(diag_pinv(F, 1, max_cond_num)));
sig_T1 = sqrt(scan_time)  * abs(sqrt(diag_pinv(F, 2, max_cond_num)));
sig_T2 = sqrt(scan_time)  * abs(sqrt(diag_pinv(F, 3, max_cond_num)));

% OLD CODE, DELETE LATER
% % Conversion to exponentials
% E1_spgr = exp(-TR_spgr ./ T1);  E2_spgr = exp(-TR_spgr ./ T2);
% E1_dess = exp(-TR_dess ./ T1);  E2_dess = exp(-TR_dess ./ T2);
% 
% % Fill in signal model gradient matrix, num_param x num_datasets
% grad_mu = []; 
% 
% % First, add columns corresponding to SPGR data
% for a = 1:length(flip_spgr)
%     grad_spgr = ...
%         [dspgr_M0s(M0s, E1_spgr, E2_spgr, flip_spgr(a), TR_spgr, TE);...
%          dspgr_T1(M0s, E1_spgr, E2_spgr, flip_spgr(a), TR_spgr, TE);...
%          dspgr_T2(M0s, E1_spgr, E2_spgr, flip_spgr(a), TR_spgr, TE)];
%     grad_mu = [grad_mu grad_spgr];
% end
% 
% % Next, add columns corresponding to DESS data
% for a = 1:length(flip_dess)
%     grad_p = ...
%         [dfp_M0s(M0s, E1_dess, E2_dess, flip_dess(a), TR_dess, TE, wf);...
%          dfp_T1(M0s, E1_dess, E2_dess, flip_dess(a), TR_dess, TE, wf);...
%          dfp_T2(M0s, E1_dess, E2_dess, flip_dess(a), TR_dess, TE, wf)];
%     grad_m = ...
%         [dfm_M0s(M0s, E1_dess, E2_dess, flip_dess(a), TR_dess, TE, wf);...
%          dfm_T1(M0s, E1_dess, E2_dess, flip_dess(a), TR_dess, TE, wf);...
%          dfm_T2(M0s, E1_dess, E2_dess, flip_dess(a), TR_dess, TE, wf)];
%     grad_mu = [grad_mu grad_p grad_m];
% end
% 
% % Construct Fisher information matrix (unit variance)
% F = grad_mu * grad_mu'; 
% 
% % Time-compensated standard deviations
% scan_time = TR_spgr * length(flip_spgr) + TR_dess * length(flip_dess);
% sig_all = sqrt(scan_time) * abs(sqrt(diag(pinv(F))));
% sig_M0s = sig_all(1);
% sig_T1 = sig_all(2);
% sig_T2 = sig_all(3);
% 
% % Alternate method of computing SDs (slower, but control cond. number)
% % max_cond_num = 1e10;
% % sig_M0s = sqrt(scan_time) * abs(sqrt(diag_pinv(F, 1, max_cond_num)));
% % sig_T1 = sqrt(scan_time) * abs(sqrt(diag_pinv(F, 2, max_cond_num)));
% % sig_T2 = sqrt(scan_time) * abs(sqrt(diag_pinv(F, 3, max_cond_num)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SPGR signal model
function sspgr = spgr(M0s, E1, E2, a, TR, TE)
sspgr = M0s .* sin(a) .* (1-E1) ./ (1-E1.*cos(a));
end

%% SPGR first derivative w.r.t. M0s
function sspgr_M0s = dspgr_M0s(M0s, E1, E2, a, TR, TE)
sspgr_M0s = sin(a) .* (1-E1) ./ (1-E1.*cos(a));
end

%% SPGR first derivative w.r.t. T1
function sspgr_T1 = dspgr_T1(M0s, E1, E2, a, TR, TE)
T1 = -TR ./ log(E1);
sspgr_E1 = M0s .* sin(a) .* (cos(a)-1) ./ (1-E1.*cos(a)).^2;
sspgr_T1 = sspgr_E1 .* E1 .* (TR ./ T1.^2);
end

%% SPGR first derivative w.r.t. T2
function sspgr_T2 = dspgr_T2(M0s, E1, E2, a, TR, TE)
T2 = -TR ./ log(E2);
sspgr_E2 = 0;
sspgr_T2 = sspgr_E2 .* E2 .* (TR ./ T2.^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DESS signal models
function sp = fp(M0s, E1, E2, a, ~, TE, wf)
sp = M0s.*tan(a.*(1.0./2.0)).*exp(pi.*5.0e-1i+TE.*wf.*1i).*...
    ((sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0).^...
    2-1.0)).*(E1-cos(a)))./(E1.*cos(a)-1.0)+1.0);
sp(isnan(sp)) = 0;
end

function sm = fm(M0s, E1, E2, a, TR, TE, wf) 
sm = -E2.^((TE.*-2.0)./TR).*M0s.*tan(a.*(1.0./2.0)).*exp(pi.*-5.0e-1i-TE.*...
    wf.*1i).*(sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./...
    (E1.*cos(a)-1.0).^2-1.0))-1.0);
sm(isnan(sm)) = 0;
end

%% DESS first derivatives w.r.t. M0s
function dsp_M0s = dfp_M0s(M0s, E1, E2, a, TR, TE, wf)
dsp_M0s = tan(a.*(1.0./2.0)).*exp(TE.*wf.*1i).*...
    ((sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0).^...
    2-1.0)).*(E1-cos(a)))./(E1.*cos(a)-1.0)+1.0).*1i;
dsp_M0s(isnan(dsp_M0s)) = 0;
end

function dsm_M0s = dfm_M0s(M0s, E1, E2, a, TR, TE, wf)
dsm_M0s = E2.^((TE.*-2.0)./TR).*tan(a.*(1.0./2.0)).*exp(TE.*wf.*-1i).*...
    (sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./...
    (E1.*cos(a)-1.0).^2-1.0))-1.0).*1i;
dsm_M0s(isnan(dsm_M0s)) = 0;
end

%% DESS first derivatives w.r.t. T1
function dsp_T1 = dfp_T1(M0s, E1, E2, a, TR, TE, wf)
T1 = -TR ./ log(E1);
dsp_E1 = M0s.*tan(a.*(1.0./2.0)).*exp(TE.*wf.*1i).*(E2.^2-1.0).*1.0./...
    sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0).^...
    2-1.0)).*(cos(a).^2-1.0).*(E1.*cos(a)-1.0).^2.*1.0./(E1.^2.*cos(a).^...
    2-E2.^2.*cos(a).^2-E1.*cos(a).*2.0-E1.^2.*E2.^2+E1.*E2.^2.*...
    cos(a).*2.0+1.0).^2.*-1i;
dsp_T1 = dsp_E1 .* E1 .* (TR ./ T1.^2);
dsp_T1(isnan(dsp_T1)) = 0;
end

function dsm_T1 = dfm_T1(M0s, E1, E2, a, TR, TE, wf)
T1 = -TR ./ log(E1);
dsm_E1 = E2.^((TE.*-2.0)./TR).*E2.^2.*M0s.*tan(a.*(1.0./2.0)).*...
    exp(TE.*wf.*-1i).*(E2.^2-1.0).*1.0./sqrt((E2.^2-1.0)./(E2.^2.*...
    (E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0).^2-1.0)).*(cos(a).^2-1.0).*...
    1.0./(E2.^2.*(E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0).^2-1.0).^2.*...
    (E1-cos(a)).*1.0./(E1.*cos(a)-1.0).^3.*-1i;
dsm_T1 = dsm_E1 .* E1 .* (TR ./ T1.^2);
dsm_T1(isnan(dsm_T1)) = 0;
end

%% DESS first derivatives w.r.t. T2
function dsp_T2 = dfp_T2(M0s, E1, E2, a, TR, TE, wf)
T2 = -TR ./ log(E2);
dsp_E2 = E2.*M0s.*tan(a.*(1.0./2.0)).*exp(TE.*wf.*1i).*(E1.^2-1.0).*...
    1.0./sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./(E1.*...
    cos(a)-1.0).^2-1.0)).*(cos(a).^2-1.0).*(E1-cos(a)).*(E1.*...
    cos(a)-1.0).*1.0./(E1.^2.*cos(a).^2-E2.^2.*cos(a).^2-E1.*...
    cos(a).*2.0-E1.^2.*E2.^2+E1.*E2.^2.*cos(a).*2.0+1.0).^2.*-1i;
dsp_T2 = dsp_E2 .* E2 .* (TR ./ T2.^2);
dsp_T2(isnan(dsp_T2)) = 0;
end

function dsm_T2 = dfm_T2(M0s, E1, E2, a, TR, TE, wf)
T2 = -TR ./ log(E2);
dsm_E2 = (E2.^(-(TE.*2.0+TR)./TR).*M0s.*TE.*tan(a.*(1.0./2.0)).*...
    exp(TE.*wf.*-1i).*(sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*...
    1.0./(E1.*cos(a)-1.0).^2-1.0))-1.0).*-2.0i)./TR-E2.^(-(TE.*...
    2.0-TR)./TR).*M0s.*tan(a.*(1.0./2.0)).*exp(TE.*wf.*-1i).*...
    (E1.^2-1.0).*1.0./sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*...
    1.0./(E1.*cos(a)-1.0).^2-1.0)).*(cos(a).^2-1.0).*(E1.*...
    cos(a)-1.0).^2.*1.0./(E1.^2.*cos(a).^2-E2.^2.*cos(a).^2-E1.*...
    cos(a).*2.0-E1.^2.*E2.^2+E1.*E2.^2.*cos(a).*2.0+1.0).^2.*1i;
dsm_T2 = dsm_E2 .* E2 .* (TR ./ T2.^2);
dsm_T2(isnan(dsm_T2)) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%