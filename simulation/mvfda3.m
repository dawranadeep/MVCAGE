%cd /storage/hpc/data/rd2nr/15.mvcage/diritchlet_cluster;
cd /Users/rdaw/Documents/MVCAGE/simulation;
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
rhs = a11^nu11 * a22^nu22/ a12^(nu11 + nu22) * gamma(nu12)/gamma(nu12 + 1) / sqrt( gamma(nu11) * gamma(nu22)  ) * sqrt( gamma(nu11 + 1) * gamma(nu22 + 1)  );
rho = unifrnd(rhs/2, 3/4 * rhs);

% Covariance function
M     = @(D, nu, a) 2^(1 - nu)/ gamma(nu) * (a * D).^nu .* besselk(nu, a * D);
cov11 = @(D) s11^2 * M(D, nu11, a11);
cov22 = @(D) s22^2 * M(D, nu22, a22);
cov12 = @(D) rho * s11 * s22 * M(D, nu12, a12);



% Sample locations
rng(100);
N  = 1000;
XX = rand(N, 1);
XX = sort(XX);


% Sample data
DD = pdist2(XX, XX);
tic;
C1  = cov11(DD); C1(1:N+ 1:end) = s11^2;
disp("First Covariance Matrix done");
toc
tic;
C2  = cov22(DD); C2(1:N+ 1:end) = s22^2;
disp("Second Covariance Matrix done");
toc
tic;
C12  = cov12(DD); C12(1:N+ 1:end) = s11 * s22 * rho;
disp("Cross Covariance Matrix done");
toc
C = [C1 C12; C12' C2];
%L = chol(C)';
nobs = 10000;
%Y = L * rand(size(L,1),nobs);
Y = mvnrnd( zeros(1, 2*N), C, nobs)';
Y1 = Y(1:N, :);
Y2 = Y((1+N):end, :);



% Plot a few simulated processes 
subplot(1,2,1); plot( XX, Y1(:,1)); xlabel("x"); ylabel("y"); title("First Sample ");
hold on;
plot( XX, Y2(:,1));
legend("First coordinate", "Second coordinate");


subplot(2,1,2); plot( XX, Y1(:,100)); xlabel("x"); ylabel("y"); title("Second Sample ");
hold on;
plot( XX, Y2(:,100));
legend("First coordinate", "Second coordinate");

sgtitle("Sample of Simulated Functional Data");

%saveas(gcf, "sample.png");
close()


%%%%%%%%%%%%%%%%    Estimated  Empirical Covariance %%%%%%%%%%%%%%%%%%%%
tic;
C_hat   = 1/(nobs-1) * Y * Y'; %cov(Y');
C1_hat  = C_hat(1:N, 1:N);
C2_hat  = C_hat((1+N):end, (1+N):end);
C12_hat = C_hat(1:N, (1+N):end);
toc





%%%%%%%%%%%%%%%%%  MVPCA %%%%%%%%%%%%%%%%%%

% Set the number of knots
n_basis = 50;

%%%%%%%%%%%%%%% First set of basis %%%%%%%%
tmp = XX * linspace(1, n_basis, n_basis);
Phi1 = [ones(N,1) sqrt(2)* sin(2* pi * tmp) sqrt(2)* cos(2* pi * tmp)];


%%%%%%%%%%%%%%%%% UKLE %%%%%%%%%%%%%%%%%%

MM = size(Phi1,2);



% UKLE1
tmp_A = zeros(MM^2, 1);
    parfor iii=1:MM^2
    i = floor((iii - 1)/MM) + 1;
    j = mod(iii - 1, MM) + 1;
    %disp([i j]);
    tmp    = C1_hat .* (Phi1(:,i) * Phi1(:,j)');
    tmp_A(iii) = trapz(XX, trapz(XX, tmp));
end
A = reshape(tmp_A, [MM MM]);
A = 1/2*(A + A');
[Evec1, Eval1]  = eigs(A, MM);      
EE1 = Phi1 * Evec1;                 % Eigenfunctions
alpha1 = zeros(nobs, MM);
parfor i=1:nobs
    if mod(i, 100) == 0
        disp(i);
    end
	alpha1(i, :) = trapz(XX, Y1(:,i) .* EE1);
end

Y1_hat = EE1 * alpha1';


%%%%%%%%%%%%% Second set of basis  %%%%%%%%%%%%


Phi2 = Phi1;


% UKLE2
tmp_A = zeros(MM^2, 1);
parfor iii=1:MM^2
    i = floor((iii - 1)/MM) + 1;
    j = mod(iii - 1, MM) + 1;
    %disp([i j]);
    tmp    = C2_hat .* (Phi2(:,i) * Phi2(:,j)');
    tmp_A(iii) = trapz(XX, trapz(XX, tmp));
end
A = reshape(tmp_A, [MM MM]);
A = 1/2*(A + A');
[Evec2, Eval2]  = eigs(A, MM);      
EE2 = Phi2 * Evec2;                 % Eigenfunctions
alpha2 = zeros(nobs, MM);
parfor i=1:nobs
    if mod(i, 100) == 0
        disp(i);
    end
alpha2(i, :) = trapz(XX, Y2(:,i) .* EE2);
end

Y2_hat = EE2 * alpha2';






% Cross KLE

K11 = Eval1;
K22 = Eval2;
K12 = zeros(MM, MM);
tmp_A = zeros(MM^2, 1);
for iii=1:MM^2
	i = floor((iii - 1)/MM) + 1;
	j = mod(iii - 1, MM) + 1;
	disp([i j]);
	tmp    = C12_hat .* (EE1(:,i) * EE2(:,j)');
	tmp_A(iii) = trapz(XX, trapz(XX, tmp));
end;
K12 = reshape(tmp_A, [MM MM]);
KK = [K11 K12; K12' K22];
[Evec, Eval] = eig(KK);



%%%%% MVKLE
EEpart1 = EE1 * Evec(1:MM, :);
EEpart2 = EE2 * Evec((MM+1):end, :);
EE = zeros(N, 2*MM, 2);
EE(:,:,1) = EEpart1;
EE(:,:,2) = EEpart2;


% scatter( XX(:,1), EE(:, 50,1),  3 );













%%%%%%%%%%%%%%%%%%%%%%%%%%% Add more pseudo points for clustering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N2 = 5000;
p_x = linspace(0, 1, N2)';

p_tmp = p_x * linspace(1, n_basis, n_basis);
p_Phi1 = [ones(N2,1) sqrt(2)* sin(2* pi * p_tmp) sqrt(2)* cos(2* pi * p_tmp)];
p_EE1    = p_Phi1 * Evec1;                 % Eigenfunctions
p_Y1     = p_EE1 * alpha1';
%plot(XX, Y1(:,7702)); hold on; plot(p_x, p_Y1(:,7702))


p_Phi2 = p_Phi1;
p_EE2    = p_Phi2 * Evec2;                 % Eigenfunctions
p_Y2     = p_EE2 * alpha2';
%plot(XX, Y2(:,7702)); hold on; plot(p_x, p_Y2(:,7702))


p_EEpart1 = p_EE1 * Evec(1:MM, :) * sqrt(Eval);
p_EEpart2 = p_EE2 * Evec((MM+1):end, :) * sqrt(Eval);
p_EE = [p_EEpart1 p_EEpart2];





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Clusters    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this_dat = [p_Y1(:,ttt)  p_Y2(:,ttt)];
this_dat = p_EE;
areas = ones(N2, 1) * 1/N2;
Z = linkage(this_dat,'ward');




rng(1);
lower = 10; %floor(N2 * 0.5/100);
upper = 200; %floor(N2 * 2/100);
cnt2 = zeros(upper - lower, 1);
cnt3 = zeros(upper - lower, 1);
parfor j= 1: (upper - lower+1)
	disp(j + lower - 1)
	%[cl, ~, sumd, ~] = kmeans(this_dat, j + lower - 1);
	cl = cluster(Z,'Maxclust',j + lower - 1);
    tmp = 0;
	for jjj = 1:max(cl)
        jj = find(cl == jjj);
        pEE = p_EE(jj, :) - mean(p_EE(jj, :));
        tmp = tmp + sum(pEE.^2, "all")/ sum(areas(jj,1));%
    end
	cnt2(j) = tmp;
	cnt3(j) = tmp/max(cl) ;
end






cnt3_dif = -cnt3(2:end) + cnt3(1:end-1);




cnt3_rel_dif = cnt3_dif./cnt3(2:end);
subplot(1,2,1); plot(lower: (upper-1), cnt3_dif);
subplot(1,2,2); plot(lower: (upper-1), cnt3_rel_dif );



best = 1+ max(find(cnt3_rel_dif > 1e-2));

rng(123);
ttt = randi([1, 10000]);

[cl, ~, sumd, ~] = kmeans([p_Y1(:,ttt)  p_Y2(:,ttt)], best + lower - 1);
yhat = zeros(N2, 2);
for j = 1:max(cl)
    jj = find(cl == j);
	yhat(jj, 1) = mean(p_Y1(jj, ttt));
	yhat(jj, 2) = mean(p_Y2(jj, ttt));
end

err = [p_Y1(:,ttt)  p_Y2(:,ttt)] - yhat;


%%%%%%%%%%%%%%%%%%%%%%%%%% Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,1);
plot(XX, Y1(:,ttt), col="blue");
hold on;
plot(p_x, yhat(:,1), col="red");  
xlabel("x");
ylabel("y");
title("First coordinate");
legend("True data", "Aggregated data");

subplot(1,2,2);  
plot(XX, Y2(:,ttt), col="blue");
hold on;
plot(p_x, yhat(:,2), col="red");
xlabel("x");
ylabel("y");
title("Second coordinate");
legend("True data", "Aggregated data");

%sgtitle("Regionalization with MVCAGE for Bivariate Matern Data");
saveas(gcf, "mvfda_3.png");
close()




%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp = 0;
cage = zeros(N2,1);
for jjj = 1:max(cl)
    jj = find(cl == jjj);
    pEE = p_EE(jj, :) - mean(p_EE(jj, :));
    cage(jj, 1) = sum(pEE.^2, "all")/ size(jj,1);
end




%subplot(1,2,1);
%plot(XX, Y1(:,ttt) - Y1_hat(:,ttt));
%hold on;
%plot(XX, Y2(:,ttt) - Y2_hat(:,ttt));
%xlabel("x");
%ylabel("y");
%title("Root-level Prediction Error");
%legend("First coordinate", "Second coordinate");


plot(lower: upper, cnt3, col = "#A2142F");
xlim([lower upper])
xlabel("Number of Regions");
ylabel("MVCAGE");
title("MVCAGE for Areal Units Number");
saveas(gcf, "cl_vs_MVCAGE.png");
close()


%subplot(1,2,2);  
plot(p_x, cage, col = "#A2142F");
xlabel("x");
ylabel("MVCAGE");
title("Aggregate-level MVCAGE");
saveas(gcf, "total_MVCAGE.png");
close()
