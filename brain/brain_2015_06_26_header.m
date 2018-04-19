%% Header Script for Brain T1/T2 Comparisons
% Requires the following files
%   1) brain_2015_06_26_ir.m
%   2) brain_2015_06_26_se.m
%   3) brain_2015_06_26_spgr_dess.m
%
% Written by: Gopal Nataraj, J. F. Nielsen, and J. A. Fessler
% Version 1.1:  2015-06-26      Based off NISTphant_2015_06_12_header.m
% Version 1.2:  2016-05-13      Adapted for manual roi selection capability

% Global parameters
fudge_factor = 1.016/0.913;     % B1 scaling factor
rs = 0;                         % ROI selection prompts
pr = 1;                         % Print .eps files
sv = 1;                         % Save data

%% T1/T2 Estimation with different SPGR/DESS combinations
% 2 SPGR, 1 DESS
n_spgr = 2;
n_dess = 1;
clear ys_im yp_im ym_im kap_reg;
brain_2015_06_26_spgr_dess();

T1_ml_2SPGR1DESS        = T1_ml;
T1_ml_mean_2SPGR1DESS   = T1_ml_mean;
T1_ml_std_2SPGR1DESS    = T1_ml_std;

T2_ml_2SPGR1DESS        = T2_ml;
T2_ml_mean_2SPGR1DESS   = T2_ml_mean;
T2_ml_std_2SPGR1DESS    = T2_ml_std;

T1_reg_2SPGR1DESS       = T1_reg;
T1_reg_mean_2SPGR1DESS  = T1_reg_mean;
T1_reg_std_2SPGR1DESS   = T1_reg_std;

T2_reg_2SPGR1DESS       = T2_reg;
T2_reg_mean_2SPGR1DESS  = T2_reg_mean;
T2_reg_std_2SPGR1DESS   = T2_reg_std;

if (sv), save('T1T2_2SPGR1DESS_brain', 'T*_2SPGR1DESS'); end

% 1 SPGR, 1 DESS
n_spgr = 1;
n_dess = 1;
clear ys_im yp_im ym_im kap_reg;
brain_2015_06_26_spgr_dess();

T1_ml_1SPGR1DESS        = T1_ml;
T1_ml_mean_1SPGR1DESS   = T1_ml_mean;
T1_ml_std_1SPGR1DESS    = T1_ml_std;

T2_ml_1SPGR1DESS        = T2_ml;
T2_ml_mean_1SPGR1DESS   = T2_ml_mean;
T2_ml_std_1SPGR1DESS    = T2_ml_std;

T1_reg_1SPGR1DESS       = T1_reg;
T1_reg_mean_1SPGR1DESS  = T1_reg_mean;
T1_reg_std_1SPGR1DESS   = T1_reg_std;

T2_reg_1SPGR1DESS       = T2_reg;
T2_reg_mean_1SPGR1DESS  = T2_reg_mean;
T2_reg_std_1SPGR1DESS   = T2_reg_std;

if (sv), save('T1T2_1SPGR1DESS_brain', 'T*_1SPGR1DESS'); end

% 0 SPGR, 2 DESS
n_spgr = 0;
n_dess = 2;
clear ys_im yp_im ym_im kap_reg;
brain_2015_06_26_spgr_dess();

T1_ml_0SPGR2DESS        = T1_ml;
T1_ml_mean_0SPGR2DESS   = T1_ml_mean;
T1_ml_std_0SPGR2DESS    = T1_ml_std;

T2_ml_0SPGR2DESS        = T2_ml;
T2_ml_mean_0SPGR2DESS   = T2_ml_mean;
T2_ml_std_0SPGR2DESS    = T2_ml_std;

T1_reg_0SPGR2DESS       = T1_reg;
T1_reg_mean_0SPGR2DESS  = T1_reg_mean;
T1_reg_std_0SPGR2DESS   = T1_reg_std;

T2_reg_0SPGR2DESS       = T2_reg;
T2_reg_mean_0SPGR2DESS  = T2_reg_mean;
T2_reg_std_0SPGR2DESS   = T2_reg_std;

if (sv), save('T1T2_0SPGR2DESS_brain', 'T*_0SPGR2DESS'); end

%% T1 Estimation from IR data
clear y_im kap_reg;
brain_2015_06_26_ir();

T1_ml_IR                = T1_ml;
T1_ml_mean_IR           = T1_ml_mean;
T1_ml_std_IR            = T1_ml_std;

T1_reg_IR               = T1_reg;
T1_reg_mean_IR          = T1_reg_mean;
T1_reg_std_IR           = T1_reg_std;

if (sv), save('T1_IR_brain', 'T*_IR'); end

%% T2 Estimation from SE data
clear y_im kap_reg;
brain_2015_06_26_se();

T2_ml_SE                = T2_ml;
T2_ml_mean_SE           = T2_ml_mean;
T2_ml_std_SE            = T2_ml_std;

T2_reg_SE               = T2_reg;
T2_reg_mean_SE          = T2_reg_mean;
T2_reg_std_SE           = T2_reg_std;

if (sv), save('T2_SE_brain', 'T*_SE'); end

%% T1/T2 Error plots for comparison
T1ideal = 0:10000;
T2ideal = 0:1000;
labels = {'WM', 'GM', 'CSF'};

% % Define SPGR/DESS ROIs (256x256)
% ctrX_ss = [109 142 101];
% ctrY_ss = [47  157 69];
% rad_ss  = [4   2    2];
% 
% % Define IR ROIs (256x256)
% ctrX_ir = [109 142 101];
% ctrY_ir = [47  157 69];
% rad_ir  = [4   2    2];
% 
% % Define SE ROIs (256x256)
% ctrX_se = [109 142 101];
% ctrY_se = [47  157 69];
% rad_se  = [4   2    2];

% Compute means and standard deviations in log domain
% Reason: this is the easiest way to make std. dev. calculations correct
logT1_reg_mean_2SPGR1DESS   = NaN(nROI, 1); 
logT1_reg_std_2SPGR1DESS    = NaN(nROI, 1);
logT1_reg_mean_1SPGR1DESS   = NaN(nROI, 1); 
logT1_reg_std_1SPGR1DESS    = NaN(nROI, 1);
logT1_reg_mean_0SPGR2DESS   = NaN(nROI, 1); 
logT1_reg_std_0SPGR2DESS    = NaN(nROI, 1);
logT1_reg_mean_IR           = NaN(nROI, 1); 
logT1_reg_std_IR            = NaN(nROI, 1);

logT2_reg_mean_2SPGR1DESS   = NaN(nROI, 1); 
logT2_reg_std_2SPGR1DESS    = NaN(nROI, 1);
logT2_reg_mean_1SPGR1DESS   = NaN(nROI, 1); 
logT2_reg_std_1SPGR1DESS    = NaN(nROI, 1);
logT2_reg_mean_0SPGR2DESS   = NaN(nROI, 1); 
logT2_reg_std_0SPGR2DESS    = NaN(nROI, 1);
logT2_reg_mean_SE           = NaN(nROI, 1); 
logT2_reg_std_SE            = NaN(nROI, 1);

for r = 1:nROI
    logT1_reg_mean_2SPGR1DESS(r)    = mean(masker(log10(T1_reg_2SPGR1DESS), roi_mask(:,:,r)));
    logT1_reg_std_2SPGR1DESS(r)     = std(masker(log10(T1_reg_2SPGR1DESS), roi_mask(:,:,r)));
    logT1_reg_mean_1SPGR1DESS(r)    = mean(masker(log10(T1_reg_1SPGR1DESS), roi_mask(:,:,r)));
    logT1_reg_std_1SPGR1DESS(r)     = std(masker(log10(T1_reg_1SPGR1DESS), roi_mask(:,:,r)));
    logT1_reg_mean_0SPGR2DESS(r)    = mean(masker(log10(T1_reg_0SPGR2DESS), roi_mask(:,:,r)));
    logT1_reg_std_0SPGR2DESS(r)     = std(masker(log10(T1_reg_0SPGR2DESS), roi_mask(:,:,r)));
    logT1_reg_mean_IR(r)            = mean(masker(log10(T1_reg_IR), roi_mask(:,:,r)));
    logT1_reg_std_IR(r)             = std(masker(log10(T1_reg_IR), roi_mask(:,:,r)));
    
    logT2_reg_mean_2SPGR1DESS(r)    = mean(masker(log10(T2_reg_2SPGR1DESS), roi_mask(:,:,r)));
    logT2_reg_std_2SPGR1DESS(r)     = std(masker(log10(T2_reg_2SPGR1DESS), roi_mask(:,:,r)));
    logT2_reg_mean_1SPGR1DESS(r)    = mean(masker(log10(T2_reg_1SPGR1DESS), roi_mask(:,:,r)));
    logT2_reg_std_1SPGR1DESS(r)     = std(masker(log10(T2_reg_1SPGR1DESS), roi_mask(:,:,r)));
    logT2_reg_mean_0SPGR2DESS(r)    = mean(masker(log10(T2_reg_0SPGR2DESS), roi_mask(:,:,r)));
    logT2_reg_std_0SPGR2DESS(r)     = std(masker(log10(T2_reg_0SPGR2DESS), roi_mask(:,:,r)));
    logT2_reg_mean_SE(r)            = mean(masker(log10(T2_reg_SE), roi_mask(:,:,r)));
    logT2_reg_std_SE(r)             = std(masker(log10(T2_reg_SE), roi_mask(:,:,r)));
    
%     [logT1_reg_mean_2SPGR1DESS(r), logT1_reg_std_2SPGR1DESS(r)] ...
%         = multiMeans(log10(T1_reg_2SPGR1DESS), [ctrX_ss(r) ctrY_ss(r)], rad_ss(r));
%     [logT1_reg_mean_1SPGR1DESS(r), logT1_reg_std_1SPGR1DESS(r)] ...
%         = multiMeans(log10(T1_reg_1SPGR1DESS), [ctrX_ss(r) ctrY_ss(r)], rad_ss(r));
%     [logT1_reg_mean_0SPGR2DESS(r), logT1_reg_std_0SPGR2DESS(r)] ...
%         = multiMeans(log10(T1_reg_0SPGR2DESS), [ctrX_ss(r) ctrY_ss(r)], rad_ss(r));
%     [logT1_reg_mean_IR(r), logT1_reg_std_IR(r)] ...
%         = multiMeans(log10(T1_reg_IR), [ctrX_ir(r) ctrY_ir(r)], rad_ir(r));
%     
%     [logT2_reg_mean_2SPGR1DESS(r), logT2_reg_std_2SPGR1DESS(r)] ...
%         = multiMeans(log10(T2_reg_2SPGR1DESS), [ctrX_ss(r) ctrY_ss(r)], rad_ss(r));
%     [logT2_reg_mean_1SPGR1DESS(r), logT2_reg_std_1SPGR1DESS(r)] ...
%         = multiMeans(log10(T2_reg_1SPGR1DESS), [ctrX_ss(r) ctrY_ss(r)], rad_ss(r));
%     [logT2_reg_mean_0SPGR2DESS(r), logT2_reg_std_0SPGR2DESS(r)] ...
%         = multiMeans(log10(T2_reg_0SPGR2DESS), [ctrX_ss(r) ctrY_ss(r)], rad_ss(r));
%     [logT2_reg_mean_SE(r), logT2_reg_std_SE(r)] ...
%         = multiMeans(log10(T2_reg_SE), [ctrX_se(r) ctrY_se(r)], rad_se(r));
end

% Create T1/T2 WM/GM ROI fill areas
T1min_roi = 800;   
T1max_roi = 1400;      
t1_xbox_roi = log10([T1min_roi T1min_roi T1max_roi T1max_roi]);
t1_ybox_roi = log10([T1min_roi T1max_roi T1max_roi T1min_roi]);

T2min_roi = 50;    
T2max_roi = 120;
t2_xbox_roi = log10([T2min_roi T2min_roi T2max_roi T2max_roi]);
t2_ybox_roi = log10([T2min_roi T2max_roi T2max_roi T2min_roi]);

% Create T1/T2 WM/GM robust-range fill areas
T1min_rob = 400;   
T1max_rob = 2000;      
t1_xbox_rob = log10([T1min_rob T1min_rob T1max_rob T1max_rob]);
t1_ybox_rob = log10([T1min_rob T1max_rob T1max_rob T1min_rob]);

T2min_rob = 40;    
T2max_rob = 200;
t2_xbox_rob = log10([T2min_rob T2min_rob T2max_rob T2max_rob]);
t2_ybox_rob = log10([T2min_rob T2max_rob T2max_rob T2min_rob]);
    
% Figure 1: T1 RLS from SPGR/DESS vs. spin-echo inversion recovery
figure; hold on;
fill(t1_xbox_rob, t1_ybox_rob, 'y');
fill(t1_xbox_roi, t1_ybox_roi, [1 0.6 0.2]);
t1_ideal = plot(log10(T1ideal), log10(T1ideal), 'k--', 'LineWidth', 1);
t1_reg_21 = errorbarxy(logT1_reg_mean_IR, logT1_reg_mean_2SPGR1DESS,...
    logT1_reg_std_IR, logT1_reg_std_2SPGR1DESS, 'Color', 'b',...
    'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 1, 'MarkerSize', 6);
t1_reg_11 = errorbarxy(logT1_reg_mean_IR, logT1_reg_mean_1SPGR1DESS,...
    logT1_reg_std_IR, logT1_reg_std_1SPGR1DESS, 'Color', 'r',...
    'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 1, 'MarkerSize', 6);
t1_reg_02 = errorbarxy(logT1_reg_mean_IR, logT1_reg_mean_0SPGR2DESS,...
    logT1_reg_std_IR, logT1_reg_std_0SPGR2DESS, 'Color', 'g',...
    'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 1, 'MarkerSize', 6);
text(logT1_reg_mean_IR, logT1_reg_mean_0SPGR2DESS, labels, 'VerticalAlignment',...
    'bottom', 'HorizontalAlignment', 'right', 'FontSize', 12);
hold off; axis([2 3.6 2 3.6]); axis square; grid on;
xlabel('SE IR log_1_0(T1) Estimates (log_1_0 ms)', 'FontSize', 14);
ylabel('SPGR/DESS log_1_0(T1) Estimates (log_1_0 ms)', 'FontSize', 14);
legend([t1_ideal, t1_reg_21(1), t1_reg_11(1), t1_reg_02(1)],...
    'Ideal', '(2,1)', '(1,1)', '(0,2)', 'Location', 'NW');
if (pr), print('-depsc', 't1_compare_brain.eps'); end;

% Figure 2: T2 RLS from SPGR/DESS vs. spin echo
figure; hold on;
fill(t2_xbox_rob, t2_ybox_rob, 'y');
fill(t2_xbox_roi, t2_ybox_roi, [1 0.6 0.2]);
t2_ideal = plot(log10(T2ideal), log10(T2ideal), 'k--', 'LineWidth', 1);
t2_reg_21 = errorbarxy(logT2_reg_mean_SE, logT2_reg_mean_2SPGR1DESS,...
    logT2_reg_std_SE, logT2_reg_std_2SPGR1DESS, 'Color', 'b',...
    'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 1, 'MarkerSize', 6);
t2_reg_11 = errorbarxy(logT2_reg_mean_SE, logT2_reg_mean_1SPGR1DESS,...
    logT2_reg_std_SE, logT2_reg_std_1SPGR1DESS, 'Color', 'r',...
    'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 1, 'MarkerSize', 6);
t2_reg_02 = errorbarxy(logT2_reg_mean_SE, logT2_reg_mean_0SPGR2DESS,...
    logT2_reg_std_SE, logT2_reg_std_0SPGR2DESS, 'Color', 'g',...
    'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 1, 'MarkerSize', 6);
text(logT2_reg_mean_SE, logT2_reg_mean_0SPGR2DESS, labels, 'VerticalAlignment',...
    'bottom', 'HorizontalAlignment', 'right', 'FontSize', 12);
hold off; axis([1 2.8 1 2.8]); axis square; grid on;
xlabel('SE log_1_0(T2) Estimates (log_1_0 ms)', 'FontSize', 14);
ylabel('SPGR/DESS log_1_0(T2) Estimates (log_1_0 ms)', 'FontSize', 14);
legend([t2_ideal, t2_reg_21(1), t2_reg_11(1), t2_reg_02(1)],...
    'Ideal', '(2,1)', '(1,1)', '(0,2)', 'Location', 'NW');
if (pr), print('-depsc', 't2_compare_brain.eps'); end;