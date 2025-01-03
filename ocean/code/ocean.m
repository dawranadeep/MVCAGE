cd /Users/rdaw/Documents/MVCAGE/ocean;
rng(22);
%parpool;
%parpool("local", 30);

YY = readmatrix("output/ocean.csv"); % First one CHA, second ROMS
Y1 = YY(:,1);
Y2 = YY(:,2);
XX = readmatrix("output/cc.csv");
%XX = YY(:, 3:4);

% Load fitted covariance matrix
C11 = readmatrix("output/C11.csv");
C12 = readmatrix("output/C12.csv");
tic;
C22 = readmatrix("output/C22.csv");
toc


N = size(YY,1);
%%%%%%%%%%%%%%%%%  MVPCA %%%%%%%%%%%%%%%%%%

% Set the number of knots in each dimension
knots = readmatrix("output/oceanKnots.csv");
%knots = rand(500, 2);
K = size(knots, 2);





%%%%%%%%%%%%%%% First set of basis %%%%%%%%
pXX = XX;
%pXX = [pXX; pXX + 0.05* randn(size(pXX))];
%pXX = [pXX; pXX + 0.05* randn(size(pXX))];

dists = pdist2(pXX, knots);
DD = pdist2(pXX, pXX);
co = 15 * min(DD + 100 *eye(4718), [], "all");
Theta1 = exp( - abs(dists)/ co );

%%%%%%%%%%%%%%%% UKLE %%%%%%%%%%%%%%%%%%

MM = size(Theta1,2);


% UKLE1
W = 1/N * Theta1' * Theta1;
pW = inv(W);
CW1 = chol(pW)';
Phi1 = exp( - abs(pdist2(XX, knots))/ co ) * CW1;               % ORthonormlalize Theta_1
W2 = 1/N * Phi1' * Phi1;
tmp = C11 * Phi1;
A2 = Phi1' * tmp;
A2 = 1/N * A2;
A2 = 1/2 * (A2 + A2');
[Evec1, Eval1]  = eigs(A2, MM);      
EE1 = Phi1 * Evec1;                 % Eigenfunctions
W3 = 1/N * EE1' * EE1;

% Check1
%W3(1:10, 1:10)
% Check2
%C1_est = 1/N * EE1 * Eval1 * EE1';
%scatter(C1_hat(:), C1_est(:), 3)


alpha1 = 1/N * Y1' * EE1;
Y1_hat = EE1 * alpha1';





%%%%%%%%%%%%% Second set of basis  %%%%%%%%%%%%
co = 5 * min(DD + 100 *eye(4718), [], "all");
Theta2 = exp( - dists.^2/2/ co^2 );
%Theta2 = [Theta2 ones(N,1)];


% UKLE2
W = 1/N * Theta2' * Theta2;
pW = inv(W);
CW2 = chol(pW)';
Phi2 = exp( - abs(pdist2(XX, knots)).^2/2/ co^2 ) * CW2;               % ORthonormlalize Theta_1
W2 = 1/N * Phi2' * Phi2;
tmp = C22 * Phi2;
A2 = Phi2' * tmp;
A2 = 1/N * A2;
A2 = 1/2 * (A2 + A2');
[Evec2, Eval2]  = eigs(A2, MM);      
EE2 = Phi2 * Evec2;                 % Eigenfunctions
W3 = 1/N * EE2' * EE2;
alpha2 = 1/N * Y2' * EE2;
Y2_hat = EE2 * alpha2';





% Cross KLE
K11 = 1/N * EE1' * C11  * EE1;
K22 = 1/N * EE2' * C22  * EE2;
K12 = 1/N * EE1' * C12 * EE2;
KK = [K11 K12; K12' K22];
%KK = [Eval1 K12; K12' Eval2];
KK = 1/2 * (KK + KK');
[Evec, Eval] = eig(KK, "vector");



%%%%% MVKLE
iii = 250:2*MM;
EEpart1 = EE1 * Evec(1:MM, iii) * diag(sqrt(Eval(iii)));
EEpart2 = EE2 * Evec((MM+1):end, iii) * diag(sqrt(Eval(iii)));
EE = [EEpart1 EEpart2];


%scatter( XX(:,1), XX(:,2),  3, EE(:, 1));






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Clusters    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[cl, cen, sumd, ddd] = kmeans(EE, 100);

%lower = 50;
%upper = 150;
%cnt = zeros(upper - lower, 1);
%parfor j= 1: (upper - lower+1)
%     disp(j + lower - 1)
%     cl = kmeans(EE, j + lower - 1);
%     sil = silhouette(EE, cl);
%     cnt(j) = mean(sil > 0.5);
%end
%plot(lower: upper, cnt);






lower = 150;
upper = 250;
cnt2 = zeros(upper - lower, 1);
cnt3 = zeros(upper - lower, 1);
parfor j= 1: (upper - lower+1)
    disp(j + lower - 1)
    [cl, ~, sumd, ~] = kmeans([Y1_hat Y2_hat], j + lower - 1);
    tmp = 0;
    for jjj = 1:max(cl)
        jj = find(cl == jjj);
        pEE = EE(jj, :) - mean(EE(jj, :));
        tmp = tmp + sum(pEE.^2, "all")/ size(jj,1);
    end
    cnt2(j) = tmp;
    cnt3(j) = tmp/max(cl) ;
end


 plot(lower: upper, cnt2);
 plot(lower: upper, cnt3);



bests = find(cnt3 == min(cnt3));
[cl, ~, sumd, ~] = kmeans([Y1_hat Y2_hat], bests + lower - 1);
yhat = zeros(N, 2);
for j = 1:max(cl)
    jj = find(cl == j);
    yhat(jj, 1) = mean(YY(jj, 1));
    yhat(jj, 2) = mean(YY(jj, 2));
end

err = YY - yhat;


bottom1 = min([YY(:,1); yhat(:,1)]);
bottom2 = min([YY(:,2); yhat(:,2)]);
top1    = max([YY(:,1); yhat(:,1)]);
top2    = max([YY(:,2); yhat(:,2)]);


% Define the first subplot
subplot(2, 2, 1);
scatter(XX(:, 1), XX(:, 2), 3, YY(:, 1), 'filled'); 
shading interp;
clim([bottom1 top1]); % Use caxis instead of clim for setting color limits
colorbar;
xlabel("Longitude");
ylabel("Latitude");
title("Original Ocean Color");

% Define the second subplot
subplot(2, 2, 2);
scatter(XX(:, 1), XX(:, 2), 3, yhat(:, 1), 'filled'); 
shading interp;
caxis([bottom1 top1]); 
colorbar;
xlabel("Longitude");
ylabel("Latitude");
title("Aggregated Ocean Color");

% Define the third subplot
subplot(2, 2, 3);
scatter(XX(:, 1), XX(:, 2), 3, YY(:, 2), 'filled'); 
shading interp;
caxis([bottom2 top2]); 
colorbar;
xlabel("Longitude");
ylabel("Latitude");
title("Original ROMS Output");

% Define the fourth subplot
subplot(2, 2, 4);
scatter(XX(:, 1), XX(:, 2), 3, yhat(:, 2), 'filled'); 
shading interp;
caxis([bottom2 top2]); 
colorbar;
xlabel("Longitude");
ylabel("Latitude");
title("Aggregated ROMS Output");
%sgtitle("MVCAGE Aggregation");
saveas(gcf, "output/ocean_new.png");
close()



