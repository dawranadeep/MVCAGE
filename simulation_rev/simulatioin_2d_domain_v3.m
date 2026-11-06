%cd /storage/hpc/data/rd2nr/15.mvcage/diritchlet_cluster;
cd /Users/rdaw/Documents/MVCAGE/simulation_rev/;

rng(22);

% =========================================================================
% 1. Generate Functional Data (Bivariate Matérn Covariance)
% =========================================================================
s11 = sqrt(2); s22 = 1;
nu11 = 0.4; nu22 = 0.5;
nu12 = 0.5 * (nu11 + nu22);
a11 = 10; a22 = 15;
a12 = 1.2 * max(a11, a22); 

rhs = (a11^nu11 * a22^nu22 / a12^(nu11 + nu22)) * ...
      (gamma(nu12) / gamma(nu12 + 1)) / ...
      sqrt(gamma(nu11) * gamma(nu22)) * sqrt(gamma(nu11 + 1) * gamma(nu22 + 1));
rho = unifrnd(rhs/2, 3/4 * rhs);

M     = @(D, nu, a) 2^(1 - nu) / gamma(nu) * (a * D).^nu .* besselk(nu, a * D);
cov11 = @(D) s11^2 * M(D, nu11, a11);
cov22 = @(D) s22^2 * M(D, nu22, a22);
cov12 = @(D) rho * s11 * s22 * M(D, nu12, a12);

% Sample 2D locations
rng(100);
N  = 1000;
XX = rand(N, 2);
DD = pdist2(XX, XX);

C1  = cov11(DD); C1(1:N+1:end) = s11^2;
C2  = cov22(DD); C2(1:N+1:end) = s22^2;
C12 = cov12(DD); C12(1:N+1:end) = s11 * s22 * rho;

C = [C1 C12; C12' C2];
nobs = 10000;
Y = mvnrnd(zeros(1, 2*N), C, nobs)';
Y1 = Y(1:N, :);
Y2 = Y((1+N):end, :);

% =========================================================================
% 2. Empirical Covariance & 2D Orthonormal Fourier Basis
% =========================================================================
C_hat   = 1/(nobs-1) * (Y * Y');
C1_hat  = C_hat(1:N, 1:N);
C2_hat  = C_hat((1+N):end, (1+N):end);
C12_hat = C_hat(1:N, (1+N):end);

n_basis = 50;
n_freq = ceil(sqrt(n_basis + 1));
[K1_grid, K2_grid] = meshgrid(0:n_freq-1, 0:n_freq-1);

K1 = K1_grid(:); K2 = K2_grid(:);
non_zero_idx = ~(K1 == 0 & K2 == 0);
K1 = K1(non_zero_idx); K1 = K1(1:n_basis);
K2 = K2(non_zero_idx); K2 = K2(1:n_basis);

arg1 = 2 * pi * (XX(:,1) * K1' + XX(:,2) * K2');
Phi1 = [ones(N,1), sqrt(2)*sin(arg1), sqrt(2)*cos(arg1)];
Phi2 = Phi1;
MM = size(Phi1, 2);

% =========================================================================
% 3. MVKLE Feature Matrix Extraction
% =========================================================================
A1 = (1/N^2) * (Phi1' * C1_hat * Phi1); A1 = 0.5 * (A1 + A1');
[Evec1, Eval1] = eigs(A1, MM);      
EE1 = Phi1 * Evec1;                 
alpha1 = (1/N) * (Y1' * EE1);      

A2 = (1/N^2) * (Phi2' * C2_hat * Phi2); A2 = 0.5 * (A2 + A2');
[Evec2, Eval2] = eigs(A2, MM);      
EE2 = Phi2 * Evec2;                 
alpha2 = (1/N) * (Y2' * EE2);      

K11 = Eval1; K22 = Eval2;
K12 = (1/N^2) * (EE1' * C12_hat * EE2);
KK = [K11 K12; K12' K22]; KK = 0.5 * (KK + KK');
[Evec, Eval] = eig(KK, 'vector');

pos_idx = Eval > 0;
Evec = Evec(:, pos_idx);
Eval = Eval(pos_idx);

% =========================================================================
% 4. Pseudo Grid & Synthetic State Partitioning (Voronoi Polygonal Domains)
% =========================================================================
N2_side = 100; 
N2 = N2_side^2; 
[p_x1, p_x2] = meshgrid(linspace(0, 1, N2_side), linspace(0, 1, N2_side));
p_XX = [p_x1(:), p_x2(:)];

p_arg = 2 * pi * (p_XX(:,1) * K1' + p_XX(:,2) * K2');
p_Phi1 = [ones(N2,1), sqrt(2)*sin(p_arg), sqrt(2)*cos(p_arg)];

p_EE1 = p_Phi1 * Evec1; p_Y1 = p_EE1 * alpha1';
p_Phi2 = p_Phi1; p_EE2 = p_Phi2 * Evec2; p_Y2 = p_EE2 * alpha2';

p_EEpart1 = p_EE1 * Evec(1:MM, :) * diag(sqrt(Eval));
p_EEpart2 = p_EE2 * Evec((MM+1):end, :) * diag(sqrt(Eval));
p_EE = [p_EEpart1 p_EEpart2];

areas = ones(N2, 1) * (1/N2);

% --- Construct Fixed "State" Polygonal Partitioning ---
rng(42);
n_states = 100; % Number of administrative states
state_seeds = rand(n_states, 2); 
[~, cl_states] = min(pdist2(p_XX, state_seeds), [], 2); % Voronoi state boundary assignments

% =========================================================================
% 5. Data-Driven Ward Hierarchical Regionalization (MVCAGE Loop)
% =========================================================================
this_dat = p_EE;
Z = linkage(this_dat, 'ward');

lower = 10;
upper = 200;
n_steps = upper - lower + 1;

cnt2 = zeros(n_steps, 1);
cnt3 = zeros(n_steps, 1);

parfor j = 1:n_steps
    k_clusters = j + lower - 1;
    cl_tmp = cluster(Z, 'Maxclust', k_clusters);
    tmp = 0;
    for jjj = 1:max(cl_tmp)
        jj = find(cl_tmp == jjj);
        pEE = p_EE(jj, :) - mean(p_EE(jj, :), 1);
        tmp = tmp + sum(pEE.^2, "all") / sum(areas(jj,1));
    end
    cnt2(j) = tmp;
    cnt3(j) = tmp / max(cl_tmp);
end

cnt3_dif = -cnt3(2:end) + cnt3(1:end-1);
cnt3_rel_dif = cnt3_dif ./ cnt3(2:end);

idx = find(cnt3_rel_dif > 1e-2);
if isempty(idx)
    best_k = lower;
else
    best_k = max(idx) + lower;
end

% Extract optimal data-driven regionalization
cl_opt = cluster(Z, 'Maxclust', best_k);

% =========================================================================
% 6. Calculate & Compare CAGE Metrics (States vs. Optimal Regionalization)
% =========================================================================

% A. Calculate CAGE for Administrative States
cage_states_map = zeros(N2, 1);
total_cage_states = 0;
for jjj = 1:max(cl_states)
    jj = find(cl_states == jjj);
    pEE = p_EE(jj, :) - mean(p_EE(jj, :), 1);
    reg_cage = sum(pEE.^2, "all") / sum(areas(jj,1));
    cage_states_map(jj) = reg_cage;
    total_cage_states = total_cage_states + reg_cage;
end
mean_cage_states = total_cage_states / max(cl_states);

% B. Calculate CAGE for Optimal Data-Driven Partition
cage_opt_map = zeros(N2, 1);
total_cage_opt = 0;
for jjj = 1:max(cl_opt)
    jj = find(cl_opt == jjj);
    pEE = p_EE(jj, :) - mean(p_EE(jj, :), 1);
    reg_cage = sum(pEE.^2, "all") / sum(areas(jj,1));
    cage_opt_map(jj) = reg_cage;
    total_cage_opt = total_cage_opt + reg_cage;
end
mean_cage_opt = total_cage_opt / max(cl_opt);

% C. Display Numerical Comparison
fprintf('\n================ REGIONALIZATION COMPARISON ================\n');
fprintf('Administrative States (Voronoi Polygons):\n');
fprintf('  - Number of Units: %d\n', max(cl_states));
fprintf('  - Total CAGE:      %.4f\n', total_cage_states);
fprintf('  - Mean CAGE:       %.4f\n\n', mean_cage_states);

fprintf('Optimal Data-Driven Regionalization:\n');
fprintf('  - Number of Units: %d\n', max(cl_opt));
fprintf('  - Total CAGE:      %.4f\n', total_cage_opt);
fprintf('  - Mean CAGE:       %.4f\n', mean_cage_opt);
fprintf('===========================================================\n');

% =========================================================================
% 7. Diagnostic Visualizations
% =========================================================================

% Plot 1: Unit Boundaries Comparison
figure('Position', [100, 100, 1000, 450]);
subplot(1,2,1);
scatter(p_XX(:,1), p_XX(:,2), 12, cl_states, 'filled');
colormap(gca, 'jet'); colorbar; xlabel("x_1"); ylabel("x_2");
title(sprintf("Administrative State Boundaries (%d Units)", max(cl_states)));

subplot(1,2,2);
scatter(p_XX(:,1), p_XX(:,2), 12, cl_opt, 'filled');
colormap(gca, 'jet'); colorbar; xlabel("x_1"); ylabel("x_2");
title(sprintf("Optimal MVCAGE Boundaries (%d Units)", max(cl_opt)));
saveas(gcf, "boundaries_comparison.png");
close();

% Plot 2: CAGE Loss Map Comparison (MAUP Evaluation)
figure('Position', [100, 100, 1000, 450]);
subplot(1,2,1);
scatter(p_XX(:,1), p_XX(:,2), 12, cage_states_map, 'filled');
colorbar; xlabel("x_1"); ylabel("x_2");
title(sprintf("State-level CAGE loss map (Mean = %.2f)", mean_cage_states));

subplot(1,2,2);
scatter(p_XX(:,1), p_XX(:,2), 12, cage_opt_map, 'filled');
colorbar; xlabel("x_1"); ylabel("x_2");
title(sprintf("Optimal CAGE Loss Map (Mean = %.2f)", mean_cage_opt));
saveas(gcf, "cage_maps_comparison.png");
close();