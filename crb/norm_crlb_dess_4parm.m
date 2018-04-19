  function [F, sig_M0s, sig_T1, sig_T2, sig_kap] ...
      = norm_crlb_dess_4parm(T1, T2, kap, wf, flip_spgr,...
      flip_dess, TR_spgr, TR_dess, TE, Sig_inv, time_comp)
% function [F, sig_M0s, sig_T1, sig_T2, sig_kap] ...
%     = norm_crlb_dess_4parm(T1, T2, kap, wf, flip_spgr,...
%     flip_dess, TR_spgr, TR_dess, TE, Sig_inv, time_comp)
%
% Inputs: 
%   T1          [1]         Spin-lattice relaxation constant
%   T2          [1]         Spin-spin relaxation constant
%   kap         [1]         Nominal flip angle compensation factor
%   wf          [1]         B0 inhomogeneity 
%   flip_spgr   [nf1]       A total nf1 of SPGR flip angle data (rad)
%   flip_dess   [nf2]       A total nf2 of DESS flip angle data (rad)
%   TR_spgr     [nf1]       SPGR repetition times (ms) 
%   TR_dess     [nf2]       DESS repetition times (ms)
%   TE          [1]         Echo time, fixed for all scans (ms)
%   Sig_inv     [M M]       Inverse of noise covariance matrix
%   time_comp   [1]         Toggle time compensation on/off
% 
% Outputs:
%   F           [P P]       Fisher info matrix (assume unit noise variance)
%   sig_M0s     [1]         Time-compensated CRLB std. dev. for M0s
%   sig_T1      [1]         Time-compensated CRLB std. dev. for T1
%   sig_T2      [1]         Time-compensated CRLB std. dev. for T2
%   sig_kap     [1]         Time-compensated CRLB std. dev. for kappa
% 
% Written by: Gopal Nataraj and Jeffrey A. Fessler; Copyright 2014
% Version 3: This CRLB is for unit M0s, allowing for SENSE normalization

% Gather constant declarations
nf1 = length(flip_spgr);
nf2 = length(flip_dess);
P = 4;                      % Number of parameters

% First, compute columns corresponding to SPGR data
grad_s = NaN(P, nf1);
for a = 1:nf1
    % Conversion to exponentials 
    E1_spgr = exp(-TR_spgr(a) ./ T1);
    E2_spgr = exp(-TR_spgr(a) ./ T2);
    
    % Append an SPGR measurement model gradient vector 
    grad_s(:,a) = ...
        [spgr(1,      E1_spgr, E2_spgr, kap, flip_spgr(a), TR_spgr(a), TE);...
         dspgr_T1(1,  E1_spgr, E2_spgr, kap, flip_spgr(a), TR_spgr(a), TE);...
         dspgr_T2(1,  E1_spgr, E2_spgr, kap, flip_spgr(a), TR_spgr(a), TE);...
         dspgr_kap(1, E1_spgr, E2_spgr, kap, flip_spgr(a), TR_spgr(a), TE)];
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
        [fp(1,      E1_dess, E2_dess, kap, flip_dess(a), TR_dess(a), TE, wf);...
         dfp_T1(1,  E1_dess, E2_dess, kap, flip_dess(a), TR_dess(a), TE, wf);...
         dfp_T2(1,  E1_dess, E2_dess, kap, flip_dess(a), TR_dess(a), TE, wf);...
         dfp_kap(1, E1_dess, E2_dess, kap, flip_dess(a), TR_dess(a), TE, wf)];
    grad_m(:,a) = ...
        [fm(1,      E1_dess, E2_dess, kap, flip_dess(a), TR_dess(a), TE, wf);...
         dfm_T1(1,  E1_dess, E2_dess, kap, flip_dess(a), TR_dess(a), TE, wf);...
         dfm_T2(1,  E1_dess, E2_dess, kap, flip_dess(a), TR_dess(a), TE, wf);...
         dfm_kap(1, E1_dess, E2_dess, kap, flip_dess(a), TR_dess(a), TE, wf)];
end

% Fill in the signal model gradient matrix, P x (M = nf1 + 2*nf2)
grad_mu = [grad_s grad_p grad_m];

% Construct Fisher information matrix (should be pure real)
F = real(grad_mu * Sig_inv * grad_mu'); 

% % Method one: standard deviations using pinv()
% if (time_comp)
%     scan_time = sum(TR_spgr) + sum(TR_dess);
%     sig_all = sqrt(scan_time) * abs(sqrt(diag(pinv(F))));
% else 
%     sig_all = abs(sqrt(diag(pinv(F))));
% end
% sig_M0s = sig_all(1);
% sig_T1 = sig_all(2);
% sig_T2 = sig_all(3);
% sig_kap = sig_all(4);

% Method two: standard deviations using condition number
max_cond_num = 1e50;
if (time_comp)
    scan_time = sum(TR_spgr) + sum(TR_dess);
    sig_M0s = sqrt(scan_time) * abs(sqrt(diag_pinv(F, 1, max_cond_num)));
    sig_T1 = sqrt(scan_time)  * abs(sqrt(diag_pinv(F, 2, max_cond_num)));
    sig_T2 = sqrt(scan_time)  * abs(sqrt(diag_pinv(F, 3, max_cond_num)));
    sig_kap = sqrt(scan_time) * abs(sqrt(diag_pinv(F, 4, max_cond_num)));
else 
    sig_M0s = abs(sqrt(diag_pinv(F, 1, max_cond_num)));
    sig_T1 = abs(sqrt(diag_pinv(F, 2, max_cond_num)));
    sig_T2 = abs(sqrt(diag_pinv(F, 3, max_cond_num)));
    sig_kap = abs(sqrt(diag_pinv(F, 4, max_cond_num)));
end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SPGR signal model
function sspgr = spgr(M0s, E1, E2, kap, a, TR, TE)
a = kap * a; 
sspgr = M0s .* sin(a) .* (1-E1) ./ (1-E1.*cos(a));
end

%% SPGR first derivative w.r.t. M0s
function sspgr_M0s = dspgr_M0s(M0s, E1, E2, kap, a, TR, TE)
a = kap * a;
sspgr_M0s = sin(a) .* (1-E1) ./ (1-E1.*cos(a));
end

%% SPGR first derivative w.r.t. T1
function sspgr_T1 = dspgr_T1(M0s, E1, E2, kap, a, TR, TE)
a = kap * a;
T1 = -TR ./ log(E1);
sspgr_E1 = M0s .* sin(a) .* (cos(a)-1) ./ (1-E1.*cos(a)).^2;
sspgr_T1 = sspgr_E1 .* E1 .* (TR ./ T1.^2);
end

%% SPGR first derivative w.r.t. T2
function sspgr_T2 = dspgr_T2(M0s, E1, E2, kap, a, TR, TE)
a = kap * a;
T2 = -TR ./ log(E2);
sspgr_E2 = 0;
sspgr_T2 = sspgr_E2 .* E2 .* (TR ./ T2.^2);
end

%% SPGR first derivative w.r.t. kap
function sspgr_kap = dspgr_kap(M0s, E1, E2, kap, a, TR, TE)
sspgr_kap = M0s.*a.*(E1-cos(a.*kap)).*1.0./(E1.*cos(a.*kap)-1.0).^2.*(E1-1.0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DESS signal models
function sp = fp(M0s, E1, E2, kap, a, TR, TE, wf)
a = kap * a;
sp = M0s.*tan(a.*(1.0./2.0)).*exp(pi.*5.0e-1i+TE.*wf.*1i).*...
    ((sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0).^...
    2-1.0)).*(E1-cos(a)))./(E1.*cos(a)-1.0)+1.0);
sp(isnan(sp)) = 0;
end

function sm = fm(M0s, E1, E2, kap, a, TR, TE, wf)
a = kap * a;
sm = -E2.^((TE.*-2.0)./TR).*M0s.*tan(a.*(1.0./2.0)).*exp(pi.*-5.0e-1i-TE.*...
    wf.*1i).*(sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./...
    (E1.*cos(a)-1.0).^2-1.0))-1.0);
sm(isnan(sm)) = 0;
end

%% DESS first derivatives w.r.t. M0s
function dsp_M0s = dfp_M0s(M0s, E1, E2, kap, a, TR, TE, wf)
a = kap * a;
dsp_M0s = tan(a.*(1.0./2.0)).*exp(TE.*wf.*1i).*...
    ((sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0).^...
    2-1.0)).*(E1-cos(a)))./(E1.*cos(a)-1.0)+1.0).*1i;
dsp_M0s(isnan(dsp_M0s)) = 0;
end

function dsm_M0s = dfm_M0s(M0s, E1, E2, kap, a, TR, TE, wf)
a = kap * a;
dsm_M0s = E2.^((TE.*-2.0)./TR).*tan(a.*(1.0./2.0)).*exp(TE.*wf.*-1i).*...
    (sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./...
    (E1.*cos(a)-1.0).^2-1.0))-1.0).*1i;
dsm_M0s(isnan(dsm_M0s)) = 0;
end

%% DESS first derivatives w.r.t. T1
function dsp_T1 = dfp_T1(M0s, E1, E2, kap, a, TR, TE, wf)
a = kap * a;
T1 = -TR ./ log(E1);
dsp_E1 = M0s.*tan(a.*(1.0./2.0)).*exp(TE.*wf.*1i).*(E2.^2-1.0).*1.0./...
    sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./(E1.*cos(a)-1.0).^...
    2-1.0)).*(cos(a).^2-1.0).*(E1.*cos(a)-1.0).^2.*1.0./(E1.^2.*cos(a).^...
    2-E2.^2.*cos(a).^2-E1.*cos(a).*2.0-E1.^2.*E2.^2+E1.*E2.^2.*...
    cos(a).*2.0+1.0).^2.*-1i;
dsp_T1 = dsp_E1 .* E1 .* (TR ./ T1.^2);
dsp_T1(isnan(dsp_T1)) = 0;
end

function dsm_T1 = dfm_T1(M0s, E1, E2, kap, a, TR, TE, wf)
a = kap * a;
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
function dsp_T2 = dfp_T2(M0s, E1, E2, kap, a, TR, TE, wf)
a = kap * a;
T2 = -TR ./ log(E2);
dsp_E2 = E2.*M0s.*tan(a.*(1.0./2.0)).*exp(TE.*wf.*1i).*(E1.^2-1.0).*...
    1.0./sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a)).^2.*1.0./(E1.*...
    cos(a)-1.0).^2-1.0)).*(cos(a).^2-1.0).*(E1-cos(a)).*(E1.*...
    cos(a)-1.0).*1.0./(E1.^2.*cos(a).^2-E2.^2.*cos(a).^2-E1.*...
    cos(a).*2.0-E1.^2.*E2.^2+E1.*E2.^2.*cos(a).*2.0+1.0).^2.*-1i;
dsp_T2 = dsp_E2 .* E2 .* (TR ./ T2.^2);
dsp_T2(isnan(dsp_T2)) = 0;
end

function dsm_T2 = dfm_T2(M0s, E1, E2, kap, a, TR, TE, wf)
a = kap * a;
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

%% DESS first derivatives w.r.t. kap
function dsp_kap = dfp_kap(M0s, E1, E2, kap, a, TR, TE, wf)
dsp_kap = M0s.*a.*exp(TE.*wf.*1i).*(tan(a.*kap.*(1.0./2.0)).^2+1.0).*...
    ((sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a.*kap)).^2.*1.0./(E1.*...
    cos(a.*kap)-1.0).^2-1.0)).*(E1-cos(a.*kap)))./(E1.*cos(a.*...
    kap)-1.0)+1.0).*5.0e-1i-M0s.*a.*exp(TE.*wf.*1i).*sin(a.*kap).*...
    tan(a.*kap.*(1.0./2.0)).*(E1.^2-1.0).*(E2.^2-1.0).*1.0./...
    sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a.*kap)).^2.*1.0./(E1.*...
    cos(a.*kap)-1.0).^2-1.0)).*(E1.*cos(a.*kap)-1.0).^2.*1.0./...
    (E1.*cos(a.*kap).*2.0+E1.^2.*E2.^2-E1.^2.*cos(a.*kap).^2+E2.^2.*...
    cos(a.*kap).^2-E1.*E2.^2.*cos(a.*kap).*2.0-1.0).^2.*1i;
dsp_kap(isnan(dsp_kap)) = 0;
end

function dsm_kap = dfm_kap(M0s, E1, E2, kap, a, TR, TE, wf)
dsm_kap = E2.^((TE.*-2.0)./TR).*M0s.*a.*exp(TE.*wf.*-1i).*(tan(a.*...
    kap.*(1.0./2.0)).^2+1.0).*(sqrt((E2.^2-1.0)./(E2.^2.*(E1-cos(a.*...
    kap)).^2.*1.0./(E1.*cos(a.*kap)-1.0).^2-1.0))-1.0).*5.0e-1i-E2.^...
    ((TE.*-2.0+TR.*2.0)./TR).*M0s.*a.*exp(TE.*wf.*-1i).*sin(a.*kap).*...
    tan(a.*kap.*(1.0./2.0)).*(E1.^2-1.0).*(E2.^2-1.0).*1.0./sqrt((E2.^...
    2-1.0)./(E2.^2.*(E1-cos(a.*kap)).^2.*1.0./(E1.*cos(a.*kap)-1.0).^...
    2-1.0)).*1.0./(E2.^2.*(E1-cos(a.*kap)).^2.*1.0./(E1.*cos(a.*...
    kap)-1.0).^2-1.0).^2.*(E1-cos(a.*kap)).*1.0./(E1.*cos(a.*kap)-1.0).^3.*1i;
dsm_kap(isnan(dsm_kap)) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
