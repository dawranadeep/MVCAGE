% cd /Users/rdaw/Documents/MVCAGE/hospital;
rng(22);

%%%%%%%%%%%%%%%%% 1. DATA SIMULATION (Replaces CSV Reading) %%%%%%%%%%%%%%%%%
N = 1000;           % Number of simulated counties
num_states = 10;    % Number of simulated states
K = 40;             % Number of knot locations for basis functions

% Simulate 2D spatial locations XX in [0,1] x [0,1]
XX = rand(N, 2);

% Assign counties to spatially contiguous "States" using k-means
state_id = kmeans(XX, num_states, 'Replicates', 5);
unique_states = unique(state_id);

% Simulate county surface areas
areas = unifrnd(0.5, 2.5, N, 1);

% Simulate correlated Bivariate Spatial Fields (Y1, Y2)
DD = pdist2(XX, XX);
s11 = 1.0; s22 = 1.0; rho = 0.65; phi = 0.25;

C11 = s11^2 * exp(-DD / phi); C11(1:N+1:end) = s11^2;
C22 = s22^2 * exp(-DD / phi); C22(1:N+1:end) = s22^2;
C12 = rho * s11 * s22 * exp(-DD / phi); C12(1:N+1:end) = rho * s11 * s22;

C = [C11, C12; C12', C22];
C = 0.5 * (C + C') + 1e-6 * eye(2*N); % Ensure symmetry and positive definiteness

Y_sim = mvnrnd(zeros(1, 2*N), C)';
Y1 = Y_sim(1:N);
Y2 = Y_sim((N+1):end);

% Standardize simulated response variables to [0, 1]
Y1 = (Y1 - min(Y1)) / (max(Y1) - min(Y1));
Y2 = (Y2 - min(Y2)) / (max(Y2) - min(Y2));

% Simulate knot locations across [0,1] x [0,1]
knots = rand(K, 2);

fprintf("Simulated Data Setup Complete:\n");
fprintf("  - Counties (N): %d\n", N);
fprintf("  - States: %d\n", num_states);
fprintf("  - Basis Knots (K): %d\n", K);

%%%%%%%%%%%%%%%%% 2. MVPCA & BASIS EXPANSION %%%%%%%%%%%%%%%%%%
pXX = XX;
pXX = [pXX; pXX + 0.05 * randn(size(pXX))];
pXX = [pXX; pXX + 0.05 * randn(size(pXX))];
dists = pdist2(pXX, knots);
co = 0.2; 
pTheta1 = exp(-abs(dists) / co);

%%%%%%%%%%%%%%%%% 3. UKLE %%%%%%%%%%%%%%%%%%
W = (1 / size(pTheta1, 1)) * (pTheta1' * pTheta1);
pW = inv(W);
CW1 = chol(pW)';
Phi = exp(-abs(pdist2(XX, knots)) / co) * CW1; 

%%%%%%%%%%%%%%% 4. GIBBS SAMPLING / G-PRIOR FITTING %%%%%%%%%%%%%
nsamp  = 50000; % Reduced iterations slightly for faster simulated runs
idx    = rand(N,1) < 1.2;
N1     = sum(idx);
g      = N;

% Posterior for Y1
PTPI   = pinv(Phi(idx, :)' * Phi(idx, :));
Vbeta  = (g / (g + 1)) * PTPI;
Ebeta  = Vbeta * Phi(idx, :)' * Y1(idx, :);
Hg     = (g / (g + 1)) * Phi(idx, :) * PTPI * Phi(idx, :)';
SSRg   = Y1(idx, :)' * (eye(sum(idx)) - Hg) * Y1(idx, :);
sig1sq = 1 ./ gamrnd(N1 / 2, 2 / SSRg, nsamp, 1);

alp1s  = zeros(nsamp, K);
parfor i = 1:nsamp
    alp1s(i, :) = mvnrnd(Ebeta, sig1sq(i) * Vbeta);
end

% Posterior for Y2
Ebeta  = Vbeta * Phi(idx, :)' * Y2(idx, :);
Hg     = (g / (g + 1)) * Phi(idx, :) * PTPI * Phi(idx, :)';
SSRg   = Y2(idx, :)' * (eye(sum(idx)) - Hg) * Y2(idx, :);
sig2sq = 1 ./ gamrnd(N1 / 2, 2 / SSRg, nsamp, 1);

alp2s  = zeros(nsamp, K);
parfor i = 1:nsamp
    alp2s(i, :) = mvnrnd(Ebeta, sig2sq(i) * Vbeta);
end

alp11s = reshape(mean(reshape(alp1s, [100, nsamp/100, K])), [nsamp/100, K]);
alp22s = reshape(mean(reshape(alp2s, [100, nsamp/100, K])), [nsamp/100, K]);

%%%%%%%%%%%%%% 5. MVKLE EIGEN-EXPANSION %%%%%%%%%%%
KK = cov([alp11s, alp22s]);
[Evec, Eval] = eig(KK, "vector");

iii = Eval > 0; % Keep positive eigenvalues
EEpart1 = Phi * Evec(1:K, iii) * diag(sqrt(Eval(iii)));
EEpart2 = Phi * Evec((K+1):end, iii) * diag(sqrt(Eval(iii)));
EE = [EEpart1, EEpart2];

%%%%%%%%%%%%%%%%%%% 6. STATE-LEVEL REGIONALIZATION & CAGE %%%%%%%%%%%%%%%%%%%%%%
state_CAGE   = zeros(num_states, 1);
county_CAGE  = zeros(N, 1);          
yhat         = zeros(N, 2);          

total_state_MVCAGE = 0;

for s = 1:num_states
    st = unique_states(s);
    jj = find(state_id == st);
    
    % Intra-state variance in coefficient space
    pEE = EE(jj, :) - mean(EE(jj, :), 1);
    
    % Compute CAGE for this state
    state_area_sum = sum(areas(jj, 1));
    state_CAGE(s)  = sum(pEE.^2, "all") / state_area_sum;
    
    county_CAGE(jj) = state_CAGE(s);
    
    % Aggregated mean spatial fields
    yhat(jj, 1) = mean(Y1(jj, :));
    yhat(jj, 2) = mean(Y2(jj, :));
    
    total_state_MVCAGE = total_state_MVCAGE + state_CAGE(s);
end

mean_state_MVCAGE = total_state_MVCAGE / num_states;

fprintf("--------------------------------------------------\n");
fprintf("Total State-Level MVCAGE Loss : %.6f\n", total_state_MVCAGE);
fprintf("Average State-Level MVCAGE Loss: %.6f\n", mean_state_MVCAGE);
fprintf("--------------------------------------------------\n");

% Display CAGE Summary Table
state_summary_tbl = table(unique_states, state_CAGE, 'VariableNames', {'Simulated_State_ID', 'CAGE'});
disp(state_summary_tbl);

%%%%%%%%%%%%%%%%%%%%%%%% 7. VISUALIZATIONS %%%%%%%%%%%%%%%%%%%%%%%%
bottom1 = min([Y1; yhat(:,1)]); top1 = max([Y1; yhat(:,1)]);
bottom2 = min([Y2; yhat(:,2)]); top2 = max([Y2; yhat(:,2)]);

figure('Position', [100 100 1000 800]);
subplot(2,2,1); scatter(XX(:,1), XX(:,2), 15, Y1(:,1), 'filled'); caxis([bottom1 top1]); colorbar; title("Simulated Y1 Field"); xlabel("x1"); ylabel("x2");
subplot(2,2,2); scatter(XX(:,1), XX(:,2), 15, yhat(:,1), 'filled'); caxis([bottom1 top1]); colorbar; title("State-Aggregated Y1 Field"); xlabel("x1"); ylabel("x2");
subplot(2,2,3); scatter(XX(:,1), XX(:,2), 15, Y2(:,1), 'filled'); caxis([bottom2 top2]); colorbar; title("Simulated Y2 Field"); xlabel("x1"); ylabel("x2");
subplot(2,2,4); scatter(XX(:,1), XX(:,2), 15, yhat(:,2), 'filled'); caxis([bottom2 top2]); colorbar; title("State-Aggregated Y2 Field"); xlabel("x1"); ylabel("x2");

figure;
scatter(XX(:,1), XX(:,2), 20, county_CAGE, 'filled'); 
colorbar; 
xlabel("x1"); ylabel("x2"); 
title("Simulated State-Level MVCAGE Map");