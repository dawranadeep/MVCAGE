%cd /storage/hpc/data/rd2nr/15.mvcage/diritchlet_cluster;
cd /Users/rdaw/Documents/MVCAGE/simulation_rev/;

rng(22);
rng(22);

% Generate functional data with a Bivariate Matern Covariance function
% Define parameters of the covariance function
s11 =  sqrt(2);
s22 = 1;
nu11 = 0.4;
nu22 = 0.5;
nu12 = 1/2*(nu11 + nu22);
a11 = 10; %0.3;
a22 = 15; %0.5;
a12 = 1.2* max(a11, a22); %0.6; has to be greater than both a11, a22
rhs = a11^nu11 * a22^nu22/ a12^(nu11 + nu22) * gamma(nu12)/gamma(nu12 + 1) / sqrt( gamma(nu11) * gamma(nu22) ) * sqrt( gamma(nu11 + 1) * gamma(nu22 + 1) );
rho = unifrnd(rhs/2, 3/4 * rhs);

% Covariance function (Isotropic Matern in 2D Euclidean Distance)
M     = @(D, nu, a) 2^(1 - nu)/ gamma(nu) * (a * D).^nu .* besselk(nu, a * D);
cov11 = @(D) s11^2 * M(D, nu11, a11);
cov22 = @(D) s22^2 * M(D, nu22, a22);
cov12 = @(D) rho * s11 * s22 * M(D, nu12, a12);

% Sample locations in 2D domain [0,1] x [0,1]
rng(100);
N  = 1000;
XX = rand(N, 2); % 2D Spatial locations (x1, x2)

% Sample data covariance matrix
DD = pdist2(XX, XX); % 2D Euclidean distance matrix
tic;
C1  = cov11(DD); C1(1:N+ 1:end) = s11^2;
disp("First Covariance Matrix done");
toc;

tic;
C2  = cov22(DD); C2(1:N+ 1:end) = s22^2;
disp("Second Covariance Matrix done");
toc;

tic;
C12  = cov12(DD); C12(1:N+ 1:end) = s11 * s22 * rho;
disp("Cross Covariance Matrix done");
toc;

C = [C1 C12; C12' C2];
nobs = 10000;
Y = mvnrnd(zeros(1, 2*N), C, nobs)';
Y1 = Y(1:N, :);
Y2 = Y((1+N):end, :);

% Plot a few simulated processes in 2D space
figure;
subplot(1,2,1); 
scatter(XX(:,1), XX(:,2), 25, Y1(:,1), 'filled'); 
colorbar; xlabel("x_1"); ylabel("x_2"); title("First Coordinate Field (Sample 1)");
subplot(1,2,2); 
scatter(XX(:,1), XX(:,2), 25, Y2(:,1), 'filled'); 
colorbar; xlabel("x_1"); ylabel("x_2"); title("Second Coordinate Field (Sample 1)");
sgtitle("Sample of Simulated 2D Spatial Functional Data");
saveas(gcf, "sample.png");
close();

%%%%%%%%%%%%%%%%    Estimated Empirical Covariance %%%%%%%%%%%%%%%%%%%%
tic;
C_hat   = 1/(nobs-1) * (Y * Y');
C1_hat  = C_hat(1:N, 1:N);
C2_hat  = C_hat((1+N):end, (1+N):end);
C12_hat = C_hat(1:N, (1+N):end);
toc;

%%%%%%%%%%%%%%%%% MVPCA %%%%%%%%%%%%%%%%%%
% Set the number of wavevector pairs
n_basis = 50;
% Generate grid of integer wavevectors (k1, k2)
n_freq = ceil(sqrt(n_basis + 1));
[K1_grid, K2_grid] = meshgrid(0:n_freq-1, 0:n_freq-1);
K1 = K1_grid(:);
K2 = K2_grid(:);
% Exclude (0,0) to prevent a zero sine column and a redundant constant cosine column
non_zero_idx = ~(K1 == 0 & K2 == 0);
K1 = K1(non_zero_idx);
K2 = K2(non_zero_idx);
% Truncate to the desired n_basis wavevectors
K1 = K1(1:n_basis);
K2 = K2(1:n_basis);
% Compute 2D Fourier arguments: 2*pi*(k1*x1 + k2*x2)
arg1 = 2 * pi * (XX(:,1) * K1' + XX(:,2) * K2');
% Build orthonormal basis matrix [1, sqrt(2)*sin(...), sqrt(2)*cos(...)]
Phi1 = [ones(N,1), sqrt(2)*sin(arg1), sqrt(2)*cos(arg1)];

%%%%%%%%%%%%%%%%% UKLE %%%%%%%%%%%%%%%%%%
MM = size(Phi1,2);
% UKLE1 (2D Integration over domain via 2D Monte Carlo / Riemann sum)
A1 = (1/N^2) * (Phi1' * C1_hat * Phi1);
A1 = 1/2*(A1 + A1');
[Evec1, Eval1]  = eigs(A1, MM);      
EE1 = Phi1 * Evec1;                 % Eigenfunctions
alpha1 = (1/N) * (Y1' * EE1);      % 2D Integral for scores
Y1_hat = EE1 * alpha1';

%%%%%%%%%%%%% Second set of basis %%%%%%%%%%%%
Phi2 = Phi1;
% UKLE2
A2 = (1/N^2) * (Phi2' * C2_hat * Phi2);
A2 = 1/2*(A2 + A2');
[Evec2, Eval2]  = eigs(A2, MM);      
EE2 = Phi2 * Evec2;                 % Eigenfunctions
alpha2 = (1/N) * (Y2' * EE2);      % 2D Integral for scores
Y2_hat = EE2 * alpha2';

% Cross KLE
K11 = Eval1;
K22 = Eval2;
K12 = (1/N^2) * (EE1' * C12_hat * EE2);
KK = [K11 K12; K12' K22];
[Evec, Eval] = eig(KK);
Eval = max(Eval, 0);

%%%%% MVKLE
EEpart1 = EE1 * Evec(1:MM, :);
EEpart2 = EE2 * Evec((MM+1):end, :);
EE = zeros(N, 2*MM, 2);
EE(:,:,1) = EEpart1;
EE(:,:,2) = EEpart2;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Add more pseudo points for clustering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N2_side = 100; 
N2 = N2_side^2; % ~4900 pseudo grid points in 2D space
[p_x1, p_x2] = meshgrid(linspace(0, 1, N2_side), linspace(0, 1, N2_side));
p_XX = [p_x1(:), p_x2(:)];
p_arg = 2 * pi * (p_XX(:,1) * K1' + p_XX(:,2) * K2');
p_Phi1 = [ones(N2,1), sqrt(2)*sin(p_arg), sqrt(2)*cos(p_arg)];
p_EE1    = p_Phi1 * Evec1;                 % Eigenfunctions
p_Y1     = p_EE1 * alpha1';
p_Phi2   = p_Phi1;
p_EE2    = p_Phi2 * Evec2;                 % Eigenfunctions
p_Y2     = p_EE2 * alpha2';
p_EEpart1 = p_EE1 * Evec(1:MM, :) * sqrt(Eval);
p_EEpart2 = p_EE2 * Evec((MM+1):end, :) * sqrt(Eval);
p_EE = [p_EEpart1 p_EEpart2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Clusters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
this_dat = p_EE;
areas = ones(N2, 1) * (1/N2);
Z = linkage(this_dat, 'ward');

rng(1);
lower = 10;
upper = 200;
cnt2 = zeros(upper - lower + 1, 1);
cnt3 = zeros(upper - lower + 1, 1);

parfor j = 1:(upper - lower + 1)
    disp(j + lower - 1);
    cl = cluster(Z, 'Maxclust', j + lower - 1);
    tmp = 0;
    for jjj = 1:max(cl)
        jj = find(cl == jjj);
        % Explicit dimension 1 prevents scalar collapse on single-point clusters
        pEE = p_EE(jj, :) - mean(p_EE(jj, :), 1);
        tmp = tmp + sum(pEE.^2, "all")/ sum(areas(jj,1));
    end
    cnt2(j) = tmp;
    cnt3(j) = tmp/max(cl);
end

cnt3_dif = -cnt3(2:end) + cnt3(1:end-1);
cnt3_rel_dif = cnt3_dif./cnt3(2:end);

figure;
subplot(1,2,1); plot(lower:(upper-1), cnt3_dif); xlabel("Number of Clusters"); title("MVCAGE Difference");
subplot(1,2,2); plot(lower:(upper-1), cnt3_rel_dif); xlabel("Number of Clusters"); title("MVCAGE Relative Difference");

% Safely extract optimal k from relative difference vector
idx = find(cnt3_rel_dif > 1e-2);
if isempty(idx)
    best_k = lower;
else
    best_k = max(idx) + lower;
end

% Extract final regionalization partition directly from Ward tree Z
cl = cluster(Z, 'Maxclust', best_k);

rng(123);
ttt = randi([1, 10000]);

yhat = zeros(N2, 2);
for j = 1:max(cl)
    jj = find(cl == j);
    yhat(jj, 1) = mean(p_Y1(jj, ttt));
    yhat(jj, 2) = mean(p_Y2(jj, ttt));
end
err = [p_Y1(:,ttt) p_Y2(:,ttt)] - yhat;

%%%%%%%%%%%%%%%%%%%%%%%%%% Spatial Data Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Spatial Data Plotting (2x2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Spatial Data Plotting (2x2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position', [100, 100, 1000, 800]);

% Determine shared color limits for Y1 and Y2 across original and aggregated fields
cLims_Y1 = [min([p_Y1(:, ttt); yhat(:,1)]), max([p_Y1(:, ttt); yhat(:,1)])];
cLims_Y2 = [min([p_Y2(:, ttt); yhat(:,2)]), max([p_Y2(:, ttt); yhat(:,2)])];

% --- Row 1, Column 1: Original First Coordinate Field ---
subplot(2,2,1);
scatter(p_XX(:,1), p_XX(:,2), 25, p_Y1(:, ttt), 'filled');
colorbar; clim(cLims_Y1); 
xlabel("x_1"); ylabel("x_2"); 
title("Original First Coordinate Field");

% --- Row 1, Column 2: Original Second Coordinate Field ---
subplot(2,2,2);
scatter(p_XX(:,1), p_XX(:,2), 25, p_Y2(:, ttt), 'filled');
colorbar; clim(cLims_Y2); 
xlabel("x_1"); ylabel("x_2"); 
title("Original Second Coordinate Field");

% --- Row 2, Column 1: Aggregated First Coordinate Field ---
subplot(2,2,3);
scatter(p_XX(:,1), p_XX(:,2), 25, yhat(:,1), 'filled');
colorbar; clim(cLims_Y1); 
xlabel("x_1"); ylabel("x_2"); 
title("Aggregated First Coordinate Field");

% --- Row 2, Column 2: Aggregated Second Coordinate Field ---
subplot(2,2,4);  
scatter(p_XX(:,1), p_XX(:,2), 25, yhat(:,2), 'filled');
colorbar; clim(cLims_Y2); 
xlabel("x_1"); ylabel("x_2"); 
title("Aggregated Second Coordinate Field");

sgtitle(sprintf("Original vs. Aggregated Spatial Fields (Realization %d)", ttt));
saveas(gcf, "mvfda_3.png");
close();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Error & Diagnostic Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cage = zeros(N2,1);
for jjj = 1:max(cl)
    jj = find(cl == jjj);
    pEE = p_EE(jj, :) - mean(p_EE(jj, :), 1);
    cage(jj, 1) = sum(pEE.^2, "all")/ size(jj,1);
end

figure;
plot(lower:upper, cnt3, 'Color', "#A2142F");
xlim([lower upper]);
xlabel("Number of Regions");
ylabel("MVCAGE");
title("MVCAGE for Areal Units Number");
saveas(gcf, "cl_vs_MVCAGE.png");
close();

figure;
scatter(p_XX(:,1), p_XX(:,2), 15, cage, 'filled');
colorbar; xlabel("x_1"); ylabel("x_2");
title("Aggregate-level MVCAGE Map");
saveas(gcf, "total_MVCAGE.png");
close();