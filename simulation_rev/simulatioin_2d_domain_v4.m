%cd /storage/hpc/data/rd2nr/15.mvcage/diritchlet_cluster;
cd /Users/rdaw/Documents/MVCAGE/simulation_rev/;

rng(22);

% =========================================================================
% 1. Generate Functional Data (Bivariate Matérn Covariance in 2D)
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

% Covariance function (Isotropic Matérn over 2D Euclidean Distance)
M     = @(D, nu, a) 2^(1 - nu) / gamma(nu) * (a * D).^nu .* besselk(nu, a * D);
cov11 = @(D) s11^2 * M(D, nu11, a11);
cov22 = @(D) s22^2 * M(D, nu22, a22);
cov12 = @(D) rho * s11 * s22 * M(D, nu12, a12);

% Sample 2D spatial locations in [0,1] x [0,1]
rng(100);
N  = 1000;
XX = rand(N, 2); 
DD = pdist2(XX, XX);

tic;
C1  = cov11(DD); C1(1:N+1:end) = s11^2;
C2  = cov22(DD); C2(1:N+1:end) = s22^2;
C12 = cov12(DD); C12(1:N+1:end) = s11 * s22 * rho;
disp("Covariance Matrices Computed.");
toc;

C = [C1 C12; C12' C2];
nobs = 10000;
Y = mvnrnd(zeros(1, 2*N), C, nobs)';
Y1 = Y(1:N, :);
Y2 = Y((1+N):end, :);

% =========================================================================
% 2. Estimated Empirical Covariance & 2D Orthonormal Fourier Basis
% =========================================================================
tic;
C_hat   = 1/(nobs-1) * (Y * Y');
C1_hat  = C_hat(1:N, 1:N);
C2_hat  = C_hat((1+N):end, (1+N):end);
C12_hat = C_hat(1:N, (1+N):end);
toc;

% 2D Basis wavevector pairs (k1, k2) excluding (0,0)
n_basis = 50;
n_freq = ceil(sqrt(n_basis + 1));
[K1_grid, K2_grid] = meshgrid(0:n_freq-1, 0:n_freq-1);

K1 = K1_grid(:); K2 = K2_grid(:);
non_zero_idx = ~(K1 == 0 & K2 == 0);
K1 = K1(non_zero_idx); K1 = K1(1:n_basis);
K2 = K2(non_zero_idx); K2 = K2(1:n_basis);

% Evaluate 2D orthonormal basis
arg1 = 2 * pi * (XX(:,1) * K1' + XX(:,2) * K2');
Phi1 = [ones(N,1), sqrt(2)*sin(arg1), sqrt(2)*cos(arg1)];
Phi2 = Phi1;
MM = size(Phi1, 2);

% =========================================================================
% 3. MVKLE Feature Matrix Extraction
% =========================================================================
% UKLE 1
A1 = (1/N^2) * (Phi1' * C1_hat * Phi1); A1 = 0.5 * (A1 + A1');
[Evec1, Eval1] = eigs(A1, MM);      
EE1 = Phi1 * Evec1;                 
alpha1 = (1/N) * (Y1' * EE1);      

% UKLE 2
A2 = (1/N^2) * (Phi2' * C2_hat * Phi2); A2 = 0.5 * (A2 + A2');
[Evec2, Eval2] = eigs(A2, MM);      
EE2 = Phi2 * Evec2;                 
alpha2 = (1/N) * (Y2' * EE2);      

% Cross KLE
K11 = Eval1; K22 = Eval2;
K12 = (1/N^2) * (EE1' * C12_hat * EE2);
KK = [K11 K12; K12' K22]; KK = 0.5 * (KK + KK');
[Evec, Eval] = eig(KK, 'vector');

pos_idx = Eval > 0;
Evec = Evec(:, pos_idx);
Eval = Eval(pos_idx);

% =========================================================================
% 4. Pseudo Grid Evaluation & MVKLE Feature Space Construction
% =========================================================================
N2_side = 100; 
N2 = N2_side^2; % 10,000 regular grid points in [0,1] x [0,1]
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

% =========================================================================
% 5. Data-Driven Hierarchical Regionalization & Optimal k Selection
% =========================================================================
this_dat = p_EE;
Z = linkage(this_dat, 'ward');

lower = 10;
upper = 200;
n_steps = upper - lower + 1; % Correct size preallocation

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

% Extract optimal data-driven partition directly from Ward tree Z
cl_opt = cluster(Z, 'Maxclust', best_k);

% =========================================================================
% 6. Construct "Administrative States" Matching best_k
% =========================================================================
% Set number of states EQUAL to the optimal cluster count derived from MVCAGE
n_states = best_k; 

rng(42);
state_seeds = rand(n_states, 2); 
[~, cl_states] = min(pdist2(p_XX, state_seeds), [], 2); % Voronoi state assignments

% =========================================================================
% 7. Calculate Aggregations & CAGE Comparisons
% =========================================================================
rng(123);
ttt = randi([1, 10000]); % Select random sample realization

% A. Aggregated Field Values for Optimal Partition
yhat = zeros(N2, 2);
for j = 1:max(cl_opt)
    jj = find(cl_opt == j);
    yhat(jj, 1) = mean(p_Y1(jj, ttt));
    yhat(jj, 2) = mean(p_Y2(jj, ttt));
end

% B. Calculate CAGE for Administrative States (k = best_k)
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

% C. Calculate CAGE for Optimal Data-Driven Partition (k = best_k)
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

% D. Print Numerical Summary
fprintf('\n================ REGIONALIZATION COMPARISON ================\n');
fprintf('Selected Optimal Units (best_k): %d\n\n', best_k);
fprintf('1. Administrative States (Voronoi Polygons):\n');
fprintf('   - Number of Units: %d\n', max(cl_states));
fprintf('   - Mean CAGE:       %.4f\n\n', mean_cage_states);
fprintf('2. Optimal Data-Driven Regionalization:\n');
fprintf('   - Number of Units: %d\n', max(cl_opt));
fprintf('   - Mean CAGE:       %.4f\n', mean_cage_opt);
fprintf('===========================================================\n');

% =========================================================================
% 8. Visualizations & Outputs
% =========================================================================

% --- Shared color limits ---
cLims_Y1 = [min([Y1(:,ttt); yhat(:,1)]), ...
            max([Y1(:,ttt); yhat(:,1)])];

cLims_Y2 = [min([Y2(:,ttt); yhat(:,2)]), ...
            max([Y2(:,ttt); yhat(:,2)])];

% --- Colormaps ---
cmap_Y1 = 'parula';
cmap_Y2 = 'cool';

% =========================================================================
% Recover pseudo-grid from p_XX
% =========================================================================
x_vals = unique(p_XX(:,1));
y_vals = unique(p_XX(:,2));

nx = length(x_vals);
ny = length(y_vals);

X_grid = reshape(p_XX(:,1), ny, nx);
Y_grid = reshape(p_XX(:,2), ny, nx);

Z_yhat1 = reshape(yhat(:,1), ny, nx);
Z_yhat2 = reshape(yhat(:,2), ny, nx);

% =========================================================================
% Figure 1: Original First Coordinate Field
% =========================================================================
figure('Position', [100, 100, 650, 550]);

ax = axes;

scatter(XX(:,1), XX(:,2), 80, Y1(:,ttt), 'filled');

clim(cLims_Y1);
colormap(ax, cmap_Y1);

axis equal tight;
xlabel("x_1", 'FontSize', 24);
ylabel("x_2", 'FontSize', 24);
title("Original first coordinate field", 'FontSize', 30);

colorbar;

saveas(gcf, "original_first_coordinate_field.png");
close();

% =========================================================================
% Figure 2: Original Second Coordinate Field
% =========================================================================
figure('Position', [100, 100, 650, 550]);

ax = axes;

scatter(XX(:,1), XX(:,2), 80, Y2(:,ttt), 'filled');

clim(cLims_Y2);
colormap(ax, cmap_Y2);

axis equal tight;
xlabel("x_1", 'FontSize', 24);
ylabel("x_2", 'FontSize', 24);
title("Original second coordinate field", 'FontSize', 30);

colorbar;

saveas(gcf, "original_second_coordinate_field.png");
close();

% =========================================================================
% Figure 3: Aggregated First Coordinate Field
% =========================================================================
figure('Position', [100, 100, 650, 550]);

ax = axes;

surf(X_grid, Y_grid, Z_yhat1, ...
    'EdgeColor', 'none');

view(2);

clim(cLims_Y1);
colormap(ax, cmap_Y1);

axis equal tight;
xlabel("x_1", 'FontSize', 24);
ylabel("x_2", 'FontSize', 24);
title("Aggregated first coordinate field", 'FontSize', 30);

colorbar;

saveas(gcf, "aggregated_first_coordinate_field.png");
close();

% =========================================================================
% Figure 4: Aggregated Second Coordinate Field
% =========================================================================
figure('Position', [100, 100, 650, 550]);

ax = axes;

surf(X_grid, Y_grid, Z_yhat2, ...
    'EdgeColor', 'none');

view(2);

clim(cLims_Y2);
colormap(ax, cmap_Y2);

axis equal tight;
xlabel("x_1", 'FontSize', 24);
ylabel("x_2", 'FontSize', 24);
title("Aggregated second coordinate field", 'FontSize', 30);

colorbar;

saveas(gcf, "aggregated_second_coordinate_field.png");
close();

% --- Plot 2: Boundary Comparison (States vs. Optimal Data-Driven) ---
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

% --- Plot 3: CAGE Loss Map Comparison (Single Shared Color Scale x 10^3) ---
figure('Position', [100, 100, 1100, 500]);

% Scale CAGE maps and mean values by 10^3 for readability
cage_states_scaled = cage_states_map / 10^3;
cage_opt_scaled    = cage_opt_map / 10^3;
mean_states_scaled = mean_cage_states / 10^3;
mean_opt_scaled    = mean_cage_opt / 10^3;

% Shared color limits computed across both scaled maps
cLims_cage = [min([cage_states_scaled; cage_opt_scaled]), ...
              max([cage_states_scaled; cage_opt_scaled])];

% Recover mesh dimensions from p_XX
x_vals = unique(p_XX(:,1));
y_vals = unique(p_XX(:,2));
nx = length(x_vals);
ny = length(y_vals);

% Reshape coordinates and values back to the mesh
X = reshape(p_XX(:,1), ny, nx);
Y = reshape(p_XX(:,2), ny, nx);

Z_states = reshape(cage_states_scaled, ny, nx);
Z_opt    = reshape(cage_opt_scaled, ny, nx);

% Tiled layout for a clean, shared colorbar setup
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Tile 1: Synthetic Regionalization ---
nexttile;
surf(X, Y, Z_states, 'EdgeColor', 'none');
view(2);
axis equal tight;
clim(cLims_cage);
xlabel("x_1", 'FontSize', 24);
ylabel("x_2", 'FontSize', 24);
title("Synthetic regionalization", 'FontSize', 30);

% --- Tile 2: Optimal Regionalization ---
nexttile;
surf(X, Y, Z_opt, 'EdgeColor', 'none');
view(2);
axis equal tight;
clim(cLims_cage);
xlabel("x_1", 'FontSize', 24);
ylabel("x_2", 'FontSize', 24);
title("Optimal regionalization", 'FontSize', 30);

% --- Single Shared Colorbar on Layout ---
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = '\times 10^3';
cb.Label.FontSize = 24;

saveas(gcf, "cage_maps_comparison.png");
close();

% --- Plot 4: MVCAGE vs. Number of Regions ---
figure;
plot(lower:upper, cnt3, 'Color', "#A2142F", 'LineWidth', 1.5);
xlim([lower upper]);
xlabel("Number of Regions");
ylabel("MVCAGE");
title("MVCAGE vs. Number of Areal Units");
saveas(gcf, "cl_vs_MVCAGE.png");
close();