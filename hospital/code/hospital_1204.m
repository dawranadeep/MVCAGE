cd /Users/rdaw/Documents/MVCAGE/hospital;
rng(22);
%parpool();

% Load the data
YY = readtable("data/county_hospital.csv"); % 2nd column median_income, 6th rating
YY = YY{:,:};
XX = YY(:, 3:4);

idx = XX(:,1) > -125;
YY = YY(idx, :);
YY(:,2) = log(YY(:,2));
Y1 = YY(:,2);
Y2 = YY(:,6);
mY1 = min(Y1); MY1 = max(Y1);
mY2 = min(Y2); MY2 = max(Y2);

% Standardize the data
Y1 = (Y1 - mY1)/ (MY1 - mY1);
Y2 = (Y2 - mY2)/ (MY2 - mY2);

XX = YY(:, 3:4);

N = size(YY,1);
%%%%%%%%%%%%%%%%%  MVPCA %%%%%%%%%%%%%%%%%%

% Read the knots dataset for basis functions
knots = readtable("data/knots.csv");
knots = knots{:,:};

m_lonlat = min(XX);
M_lonlat = max(XX);
XX = (XX - m_lonlat)./ (M_lonlat - m_lonlat);
knots = (knots - m_lonlat)./ (M_lonlat - m_lonlat);
knots = [knots; knots + 0.05* randn(size(knots))];
knots = [knots; knots + 0.05* randn(size(knots))];
K = size(knots, 1);




%%%%%%%%%%%%%%% Generate the pseudo points for GBF integration %%%%%%%%
rng(22);
pXX = XX;
pXX = [pXX; pXX + 0.05* randn(size(pXX))];
pXX = [pXX; pXX + 0.05* randn(size(pXX))];

dists = pdist2(pXX, knots);
DD = pdist2(pXX, pXX);
co = 0.2; %15 * min(DD + 100 *eye(N), [], "all");
pTheta1 = exp( - abs(dists)/ co );

%%%%%%%%%%%%%%%%% UKLE %%%%%%%%%%%%%%%%%%

MM = size(Theta1,2);


% UKLE1
W = 1/size(pTheta1,1) * pTheta1' * pTheta1;
pW = inv(W);
CW1 = chol(pW)';
Phi = exp( - abs( pdist2(XX, knots) ) / co ) * CW1;               % ORthonormlalize the actual basis function



%%%%%%%%%%%%%%% Gibbs Sampling for g prior %%%%%%%%%%%%%

alp1 = ridge(Y1, Phi, 1, 0);
YYY1 = [ones(N,1) Phi] * alp1;
alp2 = ridge(Y2, Phi, 1, 0);
YYY2 = [ones(N,1) Phi] * alp2;
sig1sq = 1/N * sum((Y1 - YYY1).^2);
sig2sq = 1/N * sum((Y2 - YYY2).^2);
subplot(1,2,1); scatter(Y1, YYY1, 3); refline([1 0])
subplot(1,2,2); scatter(Y2, YYY2, 3); refline([1 0])
alp1 = alp1(2:end);
alp2 = alp2(2:end);
sig1sq, sig2sq




%%%%%%%%%%%%%%% Integrate out sigma from posterior numerically %%%%%%%%%%%%%

nsamp  = 500000;
idx    = rand(N,1) < 1.2;
N1     = sum(idx);

g      = N;
PTPI   = pinv(Phi(idx, :)' * Phi(idx, :));
Vbeta  = (g/(g+1)) * PTPI;
Ebeta  = Vbeta * Phi(idx, :)' * Y1(idx, :);
Hg     = (g/(g+1)) * Phi(idx, :) * PTPI * Phi(idx, :)';
Yhat   = Hg*Y1(idx, :);
SSRg   = Y1(idx, :)'*  (eye(sum(idx)) - Hg)  * Y1(idx, :);
sig1sq = 1 ./gamrnd(N1/2, 2/SSRg, nsamp, 1);

alp1s  = zeros(nsamp, K);
parfor i=1:nsamp
	alp1s(i, :)  = mvnrnd(Ebeta, sig1sq(i)*Vbeta);
end






Ebeta  = Vbeta * Phi(idx, :)' * Y2(idx, :);
Hg     = (g/(g+1)) * Phi(idx, :) * PTPI * Phi(idx, :)';
Yhat   = Hg*Y2(idx, :);
SSRg   = Y2(idx, :)'*  (eye(sum(idx)) - Hg)  * Y2(idx, :);
sig2sq = 1 ./gamrnd(N1/2, 2/SSRg, nsamp, 1);

alp2s  = zeros(nsamp, K);
parfor i=1:nsamp
	alp2s(i, :)  = mvnrnd(Ebeta, sig2sq(i)*Vbeta);
end



alp11s = reshape(    mean(  reshape(alp1s, [100, nsamp/100, K])), [nsamp/100, K]   );
alp22s = reshape(    mean(  reshape(alp2s, [100, nsamp/100, K])), [nsamp/100, K]   );
alp1s  = alp11s';
alp2s  = alp22s';




%%%%%%%%%%%%%% Cross Covariance between the coefficients %%%%%%%%%%%
KK = cov([alp11s alp22s]);
[Evec, Eval] = eig(KK, "vector");

%%%%% MVKLE 
iii = Eval > 0;   % In case any negative eigenvalue
EEpart1 = Phi * Evec(1:K, iii) * diag(sqrt(Eval(iii)));
EEpart2 = Phi * Evec((K+1):end, iii) * diag(sqrt(Eval(iii)));
EE = [EEpart1 EEpart2];



%writematrix(EE, "hosp_EE.csv");









%%%%%%%%%%%%%%%%%%% Clusters %%%%%%%%%%%%%%%%%%%%%%




%EE = readmatrix("ocean_MKLE.csv"); 



this_dat = EE;
areas = YY(:,5);
Z = linkage(this_dat,'ward');    

lower = 60;
upper = 500;
cnt2 = zeros(upper - lower, 1);
cnt3 = zeros(upper - lower, 1);



parfor j= 1: (upper - lower+1)
	disp(j + lower - 1)
	%[cl, ~, sumd, ~] = kmeans(this_dat, j + lower - 1);
	cl = cluster(Z,'Maxclust',j + lower - 1);
    tmp = 0;
	for jjj = 1:max(cl)
        jj = find(cl == jjj);
        pEE = EE(jj, :) - mean(EE(jj, :));
        tmp = tmp + sum(pEE.^2, "all")/ sum(areas(jj,1));%
    end
	cnt2(j) = tmp;
	cnt3(j) = tmp/max(cl) ;
end



cnt3_dif = -cnt3(2:end) + cnt3(1:end-1);


cnt3_rel_dif = cnt3_dif./cnt3(2:end);
%subplot(1,2,1); plot(lower: (upper-1), cnt3_dif);
%subplot(1,2,2); plot(lower: (upper-1), cnt3_rel_dif );



best = 1+ max(find(cnt3_rel_dif > 1e-2));    % Find the best cluster
cl = cluster(Z,'Maxclust',best + lower - 1);   % Use this number and recluster 


yhat = zeros(N, 2);
for j = 1:max(cl)
	jj = find(cl == j);
	yhat(jj, 1) = mean(Y1(jj, :));
	yhat(jj, 2) = mean(Y2(jj, :));
end





err = ([Y1 Y2] - yhat).^2;
errhat = zeros(N, 2);
for j = 1:max(cl)
	jj = find(cl == j);
	errhat(jj, 1) = mean(err(jj, 1));
	errhat(jj, 2) = mean(err(jj, 2));
end
err = sqrt(err);
errhat = sqrt(errhat);



savedf = [YY(:,1) yhat(:,1) yhat(:,2) err(:,1) err(:,2) errhat(:,1)  errhat(:,2), cl ];
writematrix(savedf, "output/hosp.csv");


bottom1 = min([Y1  yhat(:,1)], [], 'all');
top1    = max([Y1  yhat(:,1)], [], 'all');
bottom2 = min([Y2  yhat(:,2)], [], 'all');
top2    = max([Y2  yhat(:,2)], [], 'all');

subplot(2,2,1); scatter( XX(:,1), XX(:,2), 3, Y1(:,1)  ); shading interp; caxis manual; caxis([bottom1 top1]); colorbar; xlabel("lon"); ylabel("lat"); title("Original Ocean color");
subplot(2,2,2); scatter( XX(:,1), XX(:,2), 3, yhat(:,1)); shading interp; caxis manual; caxis([bottom1 top1]); colorbar; xlabel("lon"); ylabel("lat"); title("Aggregated Ocean color");


subplot(2,2,3); scatter( XX(:,1), XX(:,2), 3, Y2(:,1)  ); shading interp; caxis manual; caxis([bottom2 top2]); colorbar; xlabel("lon"); ylabel("lat"); title("Original ROMS output");
subplot(2,2,4); scatter( XX(:,1), XX(:,2), 3, yhat(:,2)); shading interp; caxis manual; caxis([bottom2 top2]); colorbar; xlabel("lon"); ylabel("lat"); title("Aggregated ROMS output");
%sgtitle("MVCAGE Aggregation");
saveas(gcf, "output/hosp_new.png");
close()


scatter(XX(:,1), XX(:,2), 3, cl);  colorbar; xlabel("lon"); ylabel("lat"); title("MVCAGE Regions");
saveas(gcf, "output/hosp_regions.png");
close()


% Now go to the hospital_1204.R code and produce the results