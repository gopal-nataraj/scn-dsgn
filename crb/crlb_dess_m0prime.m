  function [F, sig_M0, sig_E1, sig_E2] = crlb_dess_m0prime(M0, T1, T2, ...
      wf, flip_spgr, flip_dess, TR_spgr, TR_dess, TE)
% function [F, sig_M0, sig_E1, sig_E2] = crlb_dess_m0prime(M0, T1, T2, ...
%     wf, flip_spgr, flip_dess, TR_spgr, TR_dess, TE)
% crlb_dess.m   Cramer-Rao Lower Bound (CLRB) Analysis for DESS 
% 
% Inputs: 
%   M0          [1]         Proton density, M0prime
%   T1          [1]         Spin-lattice relaxation constant
%   T2          [1]         Spin-spin relaxation constant
%   wf          [1]         B0 inhomogeneity 
%   flip_spgr   [nf1]       A total nf1 of SPGR flip angle data (rad)
%   flip_dess   [nf2]       A total nf2 of DESS flip angle data (rad)
%   TR_spgr     [1]         SPGR repetition time
%   TR_dess     [1]         DESS repetition time 
%   TE          [1]         Echo time (same for both scans)
% 
% Outputs:
%   F           [3 3]       Fisher info matrix (assume unit noise variance)
%   sig_M0      [1]         Time-compensated CRLB std. dev. for M0
%   sig_E1      [1]         Time-compensated CRLB std. dev. for E1
%   sig_E2      [1]         Time-compensated CRLB std. dev. for E2
% 
% Written by: Gopal Nataraj and Jeffrey A. Fessler; Copyright 2014

% Conversion to exponentials
E1_spgr = exp(-TR_spgr ./ T1);  E2_spgr = exp(-TR_spgr ./ T2);
E1_dess = exp(-TR_dess ./ T1);  E2_dess = exp(-TR_dess ./ T2);

% Fill in signal model gradient matrix, num_param x num_datasets
grad_mu = [];               

% First, add columns corresponding to SPGR data
for a = 1:length(flip_spgr)
    grad_spgr = ...
        [dspgr_M0(M0, E1_spgr, E2_spgr, flip_spgr(a), TR_spgr, TE);...
         dspgr_E1(M0, E1_spgr, E2_spgr, flip_spgr(a), TR_spgr, TE);...
         dspgr_E2(M0, E1_spgr, E2_spgr, flip_spgr(a), TR_spgr, TE)];
    grad_mu = [grad_mu grad_spgr];
end

% Next, add columns corresponding to DESS data
for a = 1:length(flip_dess)
    grad_p = ...
        [dfp_M0(M0, E1_dess, E2_dess, flip_dess(a), TR_dess, TE, wf);...
         dfp_E1(M0, E1_dess, E2_dess, flip_dess(a), TR_dess, TE, wf);...
         dfp_E2(M0, E1_dess, E2_dess, flip_dess(a), TR_dess, TE, wf)];
    grad_m = ...
        [dfm_M0(M0, E1_dess, E2_dess, flip_dess(a), TR_dess, TE, wf);...
         dfm_E1(M0, E1_dess, E2_dess, flip_dess(a), TR_dess, TE, wf);...
         dfm_E2(M0, E1_dess, E2_dess, flip_dess(a), TR_dess, TE, wf)];
    grad_mu = [grad_mu grad_p grad_m];
end

% Construct Fisher information matrix (unit variance)
F = grad_mu * grad_mu'; 

% Time-compensated standard deviations
scan_time = TR_spgr * length(flip_spgr) + TR_dess * length(flip_dess);
sig_all = sqrt(scan_time) * abs(sqrt(diag(pinv(F))));
sig_M0 = sig_all(1);
sig_E1 = sig_all(2);
sig_E2 = sig_all(3);

% Alternate method of computing SDs (slower, but control cond. number)
% max_cond_num = 1e10;
% sig_M0 = sqrt(scan_time) * abs(sqrt(diag_pinv(F, 1, max_cond_num)));
% sig_E1 = sqrt(scan_time) * abs(sqrt(diag_pinv(F, 2, max_cond_num)));
% sig_E2 = sqrt(scan_time) * abs(sqrt(diag_pinv(F, 3, max_cond_num)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SPGR signal model
function spgr_M0 = spgr(M0, E1, E2, a, TR, TE)
gam = TE ./ TR;
spgr_M0 = M0 .* sin(a) .* (1-E1) .* E2^gam ./ (1-E1.*cos(a));
end

%% SPGR first derivative w.r.t. M0
function sspgr_M0 = dspgr_M0(~, E1, E2, a, TR, TE)
gam = TE ./ TR;
sspgr_M0 = sin(a) .* (1-E1) .* E2^gam ./ (1-E1.*cos(a));
end

%% SPGR first derivative w.r.t. E1
function sspgr_E1 = dspgr_E1(M0, E1, E2, a, TR, TE)
gam = TE ./ TR;
sspgr_E1 = M0 .* sin(a) .* (cos(a)-1) .* E2^gam ./ (1-E1.*cos(a)).^2;
end

%% SPGR first derivative w.r.t. E2
function sspgr_E2 = dspgr_E2(M0, E1, E2, a, TR, TE)
gam = TE ./ TR;
sspgr_E2 = M0 .* sin(a) .* (1-E1) .* gam .* E2.^(gam-1) ./ (1-E1.*cos(a));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DESS signal models
function sp = fp(M0, E1, E2, a, TR, TE, wf)
sp = E2.^(TE./TR).*M0.*tan(a.*(1.0./2.0)).*exp(pi.*5.0e-1i+TE.*wf.*1i).*...
    ((sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0)...
    .^2-1.0)).*(E1-cos(a)))./(E1.*cos(a)-1.0)+1.0);
sp(isnan(sp)) = 0;
end

function sm = fm(M0, E1, E2, a, TR, TE, wf) 
sm = -E2.^(-TE./TR).*M0.*tan(a.*(1.0./2.0)).*exp(pi.*-5.0e-1i-TE.*...
    wf.*1i).*(sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./...
    (E1.*cos(a)-1.0).^2-1.0))-1.0);
sm(isnan(sm)) = 0;
end

%% DESS first derivatives w.r.t. M0
function dsp_M0 = dfp_M0(~, E1, E2, a, TR, TE, wf)
dsp_M0 = E2.^(TE./TR).*tan(a.*(1.0./2.0)).*exp(pi.*5.0e-1i+TE.*wf.*1i).*...
    ((sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0)...
    .^2-1.0)).*(E1-cos(a)))./(E1.*cos(a)-1.0)+1.0);
dsp_M0(isnan(dsp_M0)) = 0;
end

function dsm_M0 = dfm_M0(~, E1, E2, a, TR, TE, wf) 
dsm_M0 = -E2.^(-TE./TR).*tan(a.*(1.0./2.0)).*exp(pi.*-5.0e-1i-TE.*...
    wf.*1i).*(sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./...
    (E1.*cos(a)-1.0).^2-1.0))-1.0);
dsm_M0(isnan(dsm_M0)) = 0;
end

%% DESS first derivatives w.r.t. E1
function dsp_E1 = dfp_E1(M0, E1, E2, a, TR, TE, wf)
dsp_E1 = E2.^(TE./TR).*M0.*tan(a.*(1.0./2.0)).*exp(TE.*wf.*1i).*...
    (E2.^2-1.0).*1.0./sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*...
    1.0./(E1.*cos(a)-1.0).^2-1.0)).*(cos(a).^2-1.0).*(E1.*cos(a)-1.0).^...
    2.*1.0./(E1.^2.*cos(a).^2-E2.^2.*cos(a).^2-E1.*cos(a).*2.0-E1.^2.*...
    E2.^2+E1.*E2.^2.*cos(a).*2.0+1.0).^2.*-1i;
dsp_E1(isnan(dsp_E1)) = 0;
end

function dsm_E1 = dfm_E1(M0, E1, E2, a, TR, TE, wf)
dsm_E1 = E2.^(-TE./TR).*E2.^2.*M0.*tan(a.*(1.0./2.0)).*exp(TE.*wf.*-1i).*...
    (E2.^2-1.0).*1.0./sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./...
    (E1.*cos(a)-1.0).^2-1.0)).*(cos(a).^2-1.0).*1.0./(E2.^2.*...
    (E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0).^2-1.0).^2.*(E1-cos(a)).*...
    1.0./(E1.*cos(a)-1.0).^3.*-1i;
dsm_E1(isnan(dsm_E1)) = 0;
end

%% DESS first derivatives w.r.t. E2
function dsp_E2 = dfp_E2(M0, E1, E2, a, TR, TE, wf) 
dsp_E2 = (E2.^((TE-TR)./TR).*M0.*TE.*tan(a.*(1.0./2.0)).*exp(TE.*wf.*1i).*...
    ((sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0)...
    .^2-1.0)).*(E1-cos(a)))./(E1.*cos(a)-1.0)+1.0).*1i)./TR-E2.^...
    ((TE+TR)./TR).*M0.*tan(a.*(1.0./2.0)).*exp(TE.*wf.*1i).*(E1.^2-1.0)...
    .*1.0./sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./(E1.*cos(a)...
    -1.0).^2-1.0)).*(cos(a).^2-1.0).*(E1-cos(a)).*(E1.*cos(a)-1.0).*...
    1.0./(E1.^2.*cos(a).^2-E2.^2.*cos(a).^2-E1.*cos(a).*2.0-E1.^2.*...
    E2.^2+E1.*E2.^2.*cos(a).*2.0+1.0).^2.*1i;
dsp_E2(isnan(dsp_E2)) = 0;
end

function dsm_E2 = dfm_E2(M0, E1, E2, a, TR, TE, wf)
dsm_E2 = (E2.^(-(TE+TR)./TR).*M0.*TE.*tan(a.*(1.0./2.0)).*exp(TE.*wf.*-1i).*...
    (sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0).^2-...
    1.0))-1.0).*-1i)./TR-E2.^(-(TE-TR)./TR).*M0.*tan(a.*(1.0./2.0)).*...
    exp(TE.*wf.*-1i).*(E1.^2-1.0).*1.0./sqrt((E2.^2-1.0)./(E2.^2.*...
    (E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0).^2-1.0)).*(cos(a).^2-1.0).*...
    (E1.*cos(a)-1.0).^2.*1.0./(E1.^2.*cos(a).^2-E2.^2.*cos(a).^2-E1.*...
    cos(a).*2.0-E1.^2.*E2.^2+E1.*E2.^2.*cos(a).*2.0+1.0).^2.*1i;
dsm_E2(isnan(dsm_E2)) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%