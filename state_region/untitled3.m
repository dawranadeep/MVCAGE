% cd /Users/rdaw/Documents/MVCAGE/hospital;
rng(22);

fprintf("Loading real U.S. County shapefiles and geometries...\n");

%%%%%%%%%%%%%%%%% 1. REAL U.S. GEOMETRY EXTRACTION %%%%%%%%%%%%%%%%%
% Load built-in USA county map shapefile data
% Fallback: Load web shapefile if local mapping toolbox data isn't loaded
url = 'https://www2.census.gov/geo/tiger/GENZ2020/shp/cb_2020_us_county_20m.zip';
websave('counties.zip', url);
unzip('counties.zip', 'us_counties');
counties = shaperead('us_counties/cb_2020_us_county_20m.shp', 'UseGeoCoords', true);

% Filter out Non-Contiguous States/Territories (AK: 02, HI: 15, PR: 72, etc.)
exclude_state_fips = {'02', '15', '72', '60', '66', '69', '78'};
valid_idx = false(length(counties), 1);

for i = 1:length(counties)
    % Extract valid (non-NaN) boundary vertices
    lats = counties(i).Lat(~isnan(counties(i).Lat));
    lons = counties(i).Lon(~isnan(counties(i).Lon));
    
    % Skip empty geometries
    if isempty(lats) || isempty(lons)
        continue;
    end
    
    % Compute scalar centroid for the county
    mean_lat = mean(lats);
    mean_lon = mean(lons);
    
    % Extract State FIPS safely (handles both 'STATEFP' and 'StateFIPS' field names)
    if isfield(counties(i), 'STATEFP') && ~isempty(counties(i).STATEFP)
        st_fips = sprintf('%02d', str2double(counties(i).STATEFP));
    elseif isfield(counties(i), 'StateFIPS') && ~isempty(counties(i).StateFIPS)
        st_fips = sprintf('%02d', str2double(counties(i).StateFIPS));
    else
        st_fips = '';
    end
    
    % Apply logical condition with scalar values
    if ~ismember(st_fips, exclude_state_fips) && ...
       ~isnan(mean_lat) && ~isnan(mean_lon) && ...
       mean_lon > -125 % Contiguous US bounding box
        valid_idx(i) = true;
    end
end

counties = counties(valid_idx);
N = length(counties);

% Extract actual coordinates (Centroids), State FIPS, and Areas
XX_raw = zeros(N, 2);      % [Longitude, Latitude]
state_id = zeros(N, 1);    % State FIPS code
fips_code = zeros(N, 1);   % County FIPS code
areas = zeros(N, 1);       % Real county surface areas (sq km / land area)

for i = 1:N
    XX_raw(i, 1) = counties(i).Lon;
    XX_raw(i, 2) = counties(i).Lat;
    state_id(i)  = str2double(counties(i).StateFIPS);
    fips_code(i) = str2double(counties(i).FIPS);
    
    % Use land area if available, otherwise compute bounding box polygon area
    if isfield(counties(i), 'ALAND') && ~isempty(counties(i).ALAND)
        areas(i) = double(counties(i).ALAND) / 1e6; % Convert m^2 to km^2
    else
        areas(i) = (max(counties(i).Lon) - min(counties(i).Lon)) * ...
                   (max(counties(i).Lat) - min(counties(i).Lat)) * 111 * 111;
    end
end

unique_states = unique(state_id);
num_states = length(unique_states);

fprintf("Successfully loaded Real U.S. Domain:\n");
fprintf("  - Total Contiguous Counties (N): %d\n", N);
fprintf("  - Total U.S. States: %d\n", num_states);

%%%%%%%%%%%%%%%%% 2. SIMULATE BIVARIATE FIELDS OVER REAL COORDINATES %%%%%%%%%%%%%%%%%
% Normalize coordinates for stable numerical basis calculations
m_lonlat = min(XX_raw);
M_lonlat = max(XX_raw);
XX = (XX_raw - m_lonlat) ./ (M_lonlat - m_lonlat);

% Distance matrix on real geographic locations
DD = pdist2(XX, XX);
s11 = 1.0; s22 = 1.0; rho = 0.65; phi = 0.15;

C11 = s11^2 * exp(-DD / phi); C11(1:N+1:end) = s11^2;
C22 = s22^2 * exp(-DD / phi); C22(1:N+1:end) = s22^2;
C12 = rho * s11 * s22 * exp(-DD / phi); C12(1:N+1:end) = rho * s11 * s22;

C = [C11, C12; C12', C22];
C = 0.5 * (C + C') + 1e-6 * eye(2*N); % Ensure symmetry and positive definiteness

Y_sim = mvnrnd(zeros(1, 2*N), C)';
Y1 = Y_sim(1:N);
Y2 = Y_sim((N+1):end);

% Standardize response fields to [0, 1]
Y1 = (Y1 - min(Y1)) / (max(Y1) - min(Y1));
Y2 = (Y2 - min(Y2)) / (max(Y2) - min(Y2));

% Generate knot points across real normalized U.S. domain
K = 50; 
knots = XX(randperm(N, K), :); % Subsample real county centroids as knots

%%%%%%%%%%%%%%%%% 3. MVPCA & BASIS EXPANSION %%%%%%%%%%%%%%%%%%
pXX = XX;
pXX = [pXX; pXX + 0.02 * randn(size(pXX))];
pXX = [pXX; pXX + 0.02 * randn(size(pXX))];
dists = pdist2(pXX, knots);
co = 0.2; 
pTheta1 = exp(-abs(dists) / co);

W = (1 / size(pTheta1, 1)) * (pTheta1' * pTheta1);
pW = inv(W);
CW1 = chol(pW)';
Phi = exp(-abs(pdist2(XX, knots)) / co) * CW1; 

%%%%%%%%%%%%%%%%% 4. GIBBS SAMPLING / G-PRIOR FITTING %%%%%%%%%%%%%%%%%
nsamp  = 50000;
idx    = rand(N,1) < 1.2;
N1     = sum(idx);
g      = N;

% Posterior Y1
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

% Posterior Y2
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

%%%%%%%%%%%%%%%%% 5. MVKLE EIGEN-EXPANSION %%%%%%%%%%%%%%%%%
KK = cov([alp11s, alp22s]);
[Evec, Eval] = eig(KK, "vector");

iii = Eval > 0;
EEpart1 = Phi * Evec(1:K, iii) * diag(sqrt(Eval(iii)));
EEpart2 = Phi * Evec((K+1):end, iii) * diag(sqrt(Eval(iii)));
EE = [EEpart1, EEpart2];

%%%%%%%%%%%%%%%%% 6. REAL STATE-LEVEL REGIONALIZATION & CAGE %%%%%%%%%%%%%%%%%
state_CAGE   = zeros(num_states, 1);
county_CAGE  = zeros(N, 1);          
yhat         = zeros(N, 2);          

total_state_MVCAGE = 0;

for s = 1:num_states
    st = unique_states(s);
    jj = find(state_id == st);
    
    % Intra-state variance in coefficient space
    pEE = EE(jj, :) - mean(EE(jj, :), 1);
    
    % CAGE for actual state
    state_area_sum = sum(areas(jj, 1));
    state_CAGE(s)  = sum(pEE.^2, "all") / state_area_sum;
    
    county_CAGE(jj) = state_CAGE(s);
    
    % Aggregated mean fields across actual state boundary
    yhat(jj, 1) = mean(Y1(jj, :));
    yhat(jj, 2) = mean(Y2(jj, :));
    
    total_state_MVCAGE = total_state_MVCAGE + state_CAGE(s);
end

mean_state_MVCAGE = total_state_MVCAGE / num_states;

fprintf("--------------------------------------------------\n");
fprintf("Total U.S. State-Level MVCAGE Loss : %.6f\n", total_state_MVCAGE);
fprintf("Average State-Level MVCAGE Loss    : %.6f\n", mean_state_MVCAGE);
fprintf("--------------------------------------------------\n");

% Summary Table
state_summary_tbl = table(unique_states, state_CAGE, 'VariableNames', {'State_FIPS', 'CAGE'});
disp(head(state_summary_tbl, 10));

%%%%%%%%%%%%%%%%% 7. GEOGRAPHIC MAP VISUALIZATIONS %%%%%%%%%%%%%%%%%
bottom1 = min([Y1; yhat(:,1)]); top1 = max([Y1; yhat(:,1)]);
bottom2 = min([Y2; yhat(:,2)]); top2 = max([Y2; yhat(:,2)]);

% Map 1: Original vs Aggregated Fields on US Map
figure('Position', [100 100 1100 800]);
subplot(2,2,1); scatter(XX_raw(:,1), XX_raw(:,2), 6, Y1(:,1), 'filled'); caxis([bottom1 top1]); colorbar; title("Simulated County Field Y1"); xlabel("Longitude"); ylabel("Latitude");
subplot(2,2,2); scatter(XX_raw(:,1), XX_raw(:,2), 6, yhat(:,1), 'filled'); caxis([bottom1 top1]); colorbar; title("U.S. State-Aggregated Y1"); xlabel("Longitude"); ylabel("Latitude");
subplot(2,2,3); scatter(XX_raw(:,1), XX_raw(:,2), 6, Y2(:,1), 'filled'); caxis([bottom2 top2]); colorbar; title("Simulated County Field Y2"); xlabel("Longitude"); ylabel("Latitude");
subplot(2,2,4); scatter(XX_raw(:,1), XX_raw(:,2), 6, yhat(:,2), 'filled'); caxis([bottom2 top2]); colorbar; title("U.S. State-Aggregated Y2"); xlabel("Longitude"); ylabel("Latitude");

% Map 2: Actual U.S. State CAGE Loss Distribution Map
figure('Position', [150 150 900 600]);
scatter(XX_raw(:,1), XX_raw(:,2), 10, county_CAGE, 'filled'); 
colorbar; 
xlabel("Longitude"); ylabel("Latitude"); 
title("Calculated CAGE Loss When Aggregating Counties to U.S. States");