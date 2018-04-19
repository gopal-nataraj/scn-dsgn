
flip1 = 25:10:75;
flip2 = 40:10:80;
rmse_gm = NaN(length(flip1), length(flip2));
rmse_wm = NaN(length(flip1), length(flip2));

for ii = 1:length(flip1)
    for jj = 1:length(flip2);
        flip = [flip1(ii) flip2(jj) 85]' * pi/180;
        nf = length(flip);
        TR = 20;                                % ms
        TE = 5;                                 % ms
        wf_true = 0;                            % Off-resonance
        SNR = 40;                               % dB

        % Forward model: make true data, using M0s_true
        yp_true = NaN(nx, ny, nf);
        ym_true = NaN(nx, ny, nf);
        for a = 1:nf
            v1 = (1 - E1_true * cos(flip(a))) ./ (E1_true - cos(flip(a)));
            [S1true, S2true] = dess_fun_M0star(TR, TE, M0s_true, wf_true, flip(a), v1, T2_true);
            yp_true(:,:,a) = fft2(S1true);
            ym_true(:,:,a) = fft2(S2true);
        end

        % Add complex white gaussian noise
        sigma_p = exp(-SNR/20) * norm(yp_true(:)) / sqrt(2*numel(yp_true));
        sigma_m = exp(-SNR/20) * norm(ym_true(:)) / sqrt(2*numel(ym_true));
        yp = yp_true + sigma_p * (randn(size(yp_true)) + 1i * randn(size(yp_true)));
        ym = ym_true + sigma_m * (randn(size(ym_true)) + 1i * randn(size(ym_true)));
        printm('snr_p = %g dB', 20*log(norm(yp_true(:)) / norm(col(yp_true-yp))));
        printm('snr_m = %g dB', 20*log(norm(ym_true(:)) / norm(col(ym_true-ym))));

        % Take noisy data back to image domain
        % Since simulating only body coil, sqrt(SOS) is just the magnitude
        yp_im = ifft2(yp);
        ym_im = ifft2(ym);

        % Method-of-Moments T2 Estimation
        if (nf == 1)
            T2_mom = -2*(TR-TE) ./ log(abs(ym_im ./ yp_im));
        else
            E2squared = NaN(nx, ny);
            for xx = 1:nx
                for yy = 1:ny 
                    Y = squeeze(ym_im(xx,yy,:));
                    A = squeeze(yp_im(xx,yy,:));
                    C = [A Y];
                    sig_small = min(svds(C'*C));
                    E2squared(xx,yy) = abs((A'*A - sig_small.^2) \ (A' * Y));
                end
            end
            T2_mom = -2*(TR-TE) ./ log(E2squared);
        end
        
        rmse_gm(ii,jj) = rmse(T2_mom(mask_gm), T2_true(mask_gm));
        rmse_wm(ii,jj) = rmse(T2_mom(mask_wm), T2_true(mask_wm));
    end
end

figure; im(flip1, flip2, rmse_gm);
figure; im(flip1, flip2, rmse_wm);
figure; im(flip1, flip2, rmse_gm + rmse_wm);
        
        