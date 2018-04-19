%% Test script to understand/approximate noncentral chi distribution
%% Create empirical distribution in the real case
nc = 8;
ns = 100000;

mu = linspace(0.1, 0.8, nc);            % Coil-dependent means
sd = 1e-3 * ones(1, nc);                % Same standard deviation 
X = NaN(ns, nc);
for i = 1:nc
    X(:,i) = mu(i) + sd(i) .* randn(ns,1);
end

Z = sqrt(sum((X ./ repmat(sd, [ns 1])).^2, 2));
figure; hist(Z, 50); 
mean_emp = mean(Z);
var_emp = var(Z);

% Compare with theoretical distribution
lam = sqrt(sum((mu ./ sd).^2)); 
mean_theo = sqrt(pi/2) * mfun('L', 1/2, nc/2-1, -lam.^2/2); 
var_theo = nc + lam.^2 - mean_theo.^2;

%% Create empirical distribution in the complex case 
nc = 8;
ns = 100000;

mu = linspace(0, 0, nc);
%mu = linspace(0.1, 0.8, nc) .* exp(1i * linspace(0.1, 0.8, nc));
sd_r = (1e-3/sqrt(2)) * ones(1, nc);
sd_i = sd_r;                            % Circularly symmetric assumption
X = NaN(ns, nc);
for i = 1:nc
    X(:,i) = mu(i) + (sd_r(i).*randn(ns,1)) + (1i.*sd_i(i).*randn(ns,1));
end

Z = sqrt(sum(abs(X ./ repmat(sd_r, [ns 1])).^2, 2));
figure; hist(Z, 50);
mean_emp = mean(Z); 
var_emp = var(Z); 

% Compare with theoretical distribution 
lam = sqrt(sum((abs(mu) ./ sd_r).^2));
mean_theo = sqrt(pi/2) * mfun('L', 1/2, nc-1, -lam.^2/2);
var_theo = 2*nc + lam.^2 - mean_theo.^2;