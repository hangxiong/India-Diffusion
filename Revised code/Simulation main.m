%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main_models for experience effect and externality effect

% Outline

% 1. Select villages to consider
% 2. Select moments
% 3. Select time vector and number of repetitions per trial
% 4. Select parameter grid
% 5. Load data
% 6. Logistic fit to get coefficients for covariates and the constant
% 7. RUNNING THE MODEL
% 8. RUNNING THE AGGREGATOR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
% Adjusting so this runs off dropbox
cd ..
clear all

location = pwd;
addpath(genpath(location));

%% Parameters

% 0. Replicate the original study?
replicateOriginal = 1;

% 1. Select villages to consider
vills = [1:4,6,9, 12, 15, 19:21, 23:25, 29, 31:33, 36, 39, 42, 43, 45:48, 50:52, 55, 57, 59:60, 62, 64:65, 67:68, 70:73, 75];
G = length(vills); % Number of graphs

% 2. Select moments
version = 1;

switch version
    case 1
        m = 5;
    case 2
        m = 3;
    case 3
        m = 3;
    case 4
        m = 3;
end

% 3. Select time vector and number of repetitions per trial
S = 75; % Number of simulations
timeVector = 'trimesters';

TMonths = [31 35 15 35 13 2 32 5 31 35 31 29 19 22 25 25 23 23 24 25 26 24 17 16 17 13 19 20 20 19 22 14 12 15 10 19 18 18 19 19 19 17 17]; % Months

switch timeVector
    case 'months'
        T = TMonths + 1;
    case 'quarters'
        T = ceil(TMonths./3) + 1;  %Quarters have 3 months in them
    case 'trimesters'
        T = ceil(TMonths./4) + 1; %Trimesters have 4 months in them
end
%T = T(1:5);
assert(G == numel(T))


% 4. Select parameter grid
%lambda = [(-1 : 0.1 : -0.3), (-0.25 : 0.05 : 0.3), (0.4 : 0.1 : 1)];
%qN = [(0 : 0.05 : 0.5), (0.6 : 0.1: 1)];
%qP = [(0 : 0.05 : 0.5), (0.6 : 0.1: 1)];
%threshold = [(0.5 : 0.05 : 1)];

lambda = [0.1 0.2];
qN = [0.02 0.04];
qP = [0.05 0.1];
threshold = [0.6 0.7];

% 5. Load data
relative = 1;

%% Pre-allocation
inGiant = cell(G,1);
leaders = cell(G,1);
TakeUp = cell(G,1);
TakingLeaders = cell(G,1);
ZLeaders = cell(G,1);
Covars = [];
Outcome = [];
EmpRate = zeros(G,1);
rdist = cell(G,1);
dist = cell(G,1);
hermits = cell(G,1);
W = cell(G,1);
Z = cell(G,1);
EOmicronX = cell(G, 1);
EOmicronY = cell(G, 1);
EOmicronM = cell(G, 1);
DOmicronX = cell(G, 1);
DOmicronY = cell(G, 1);
DOmicronM = cell(G, 1);
NOmicronX = cell(G, 1);
NOmicronY = cell(G, 1);
NOmicronM = cell(G, 1);

% Load the household connection adjacency matrix.
% networks of all relationships
X = load(['network data/AdjacencyMatrics_All.mat']);
X = X.adjX;
% networks of non borrowing and lending relationships
Y = load(['network data/AdjacencyMatrics_NonMoney.mat']);
Y = Y.adjY;
% networks of borrowing and lending relationships
M = load(['network data/AdjacencyMatrics_Money.mat']);
M = M.adjM;

if(replicateOriginal)
    Y = X;
    M = X;
end

%% Construct data
counter = 0;
for vilnum = vills
    counter = counter + 1;
    
    N = length(X{counter});
    
    if relative==0
        EOmicronX{counter} = dlmread(['./network data/EOmicron_abs_X' num2str(vilnum) '.csv']);
        EOmicronY{counter} = dlmread(['./network data/EOmicron_abs_Y' num2str(vilnum) '.csv']);
        EOmicronM{counter} = dlmread(['./network data/EOmicron_abs_M' num2str(vilnum) '.csv']);
        DOmicronX{counter} = dlmread(['./network data/DOmicron_abs_X' num2str(vilnum) '.csv']);
        DOmicronY{counter} = dlmread(['./network data/DOmicron_abs_Y' num2str(vilnum) '.csv']);
        DOmicronM{counter} = dlmread(['./network data/DOmicron_abs_M' num2str(vilnum) '.csv']);
        NOmicronX{counter} = dlmread(['./network data/NOmicron_abs_X' num2str(vilnum) '.csv']);
        NOmicronY{counter} = dlmread(['./network data/NOmicron_abs_Y' num2str(vilnum) '.csv']);
        NOmicronM{counter} = dlmread(['./network data/NOmicron_abs_M' num2str(vilnum) '.csv']);
    elseif relative==1,
        EOmicronX{counter} = dlmread(['./network data/EOmicron_rel_X' num2str(vilnum) '.csv']);
        EOmicronY{counter} = dlmread(['./network data/EOmicron_rel_Y' num2str(vilnum) '.csv']);
        EOmicronM{counter} = dlmread(['./network data/EOmicron_rel_M' num2str(vilnum) '.csv']);
        DOmicronX{counter} = dlmread(['./network data/DOmicron_rel_X' num2str(vilnum) '.csv']);
        DOmicronY{counter} = dlmread(['./network data/DOmicron_rel_Y' num2str(vilnum) '.csv']);
        DOmicronM{counter} = dlmread(['./network data/DOmicron_rel_M' num2str(vilnum) '.csv']);
        NOmicronX{counter} = dlmread(['./network data/NOmicron_rel_X' num2str(vilnum) '.csv']);
        NOmicronY{counter} = dlmread(['./network data/NOmicron_rel_Y' num2str(vilnum) '.csv']);
        NOmicronM{counter} = dlmread(['./network data/NOmicron_rel_M' num2str(vilnum) '.csv']);
    end
    
    if(replicateOriginal)
        EOmicronY = EOmicronX;
        EOmicronM = EOmicronX;
        DOmicronY = DOmicronX;
        DOmicronM = DOmicronX;
        NOmicronY = NOmicronX;
        NOmicronM = NOmicronX;
    end
    
    % Load the Leader data
    templeaders = load(['./demographic data/HHhasALeader' num2str(vilnum) '.csv']);
    leaders{counter} = templeaders(:,2);
    
    % Load the Take-Up data
    TakeUp{counter} = load(['./demographic data/MF' num2str(vilnum) '.csv']);
    EmpRate(counter) = mean(TakeUp{counter}(~leaders{counter})); % leaders are excluded
    
    % Load the giant component data
    inGiant{counter} = load(['./demographic data/inGiant' num2str(vilnum) '.csv']);
    
    % Generate hermits
    d = sum(X{counter},2); % Sum up along columns
    hermits{counter}=(d==0);
    
    % Load the Covariates data
    W{counter} = load(['./demographic data/hhcovariates' num2str(vilnum) '.csv']);
    
    % Which covariates to use - for instance we want to add a PCA
    Z{counter} = [W{counter}(:,1:6)]; % for instance, take the first covariate only
    
    % Prune other stats, i.e., drop off notes that are not in giants, etc.
    leaders{counter} = leaders{counter}(logical(inGiant{counter}));
    TakeUp{counter} = TakeUp{counter}(logical(inGiant{counter}));
    Z{counter} = Z{counter}(logical(inGiant{counter}),:);
    
    TakingLeaders{counter} = TakeUp{counter}(logical(leaders{counter}));
    ZLeaders{counter} = Z{counter}(logical(leaders{counter}),:);
    Outcome = [Outcome; TakingLeaders{counter}];
    Covars = [Covars; ZLeaders{counter}];
    
    % Second neighbors
    Sec{counter} = (X{counter}^2>0);
    for i=1:length(X{counter})
        Sec{counter}(i,i) = 0;
    end
    Sec{counter}=(Sec{counter}-X{counter} > 0);
end


% 6. Logistic fit to get coefficients for covariates and the constant
[Betas, dev, stats] = glmfit(Covars,Outcome,'binomial','link','logit');
[Betas'; stats.se'; stats.p'];

tic;

% 7. RUNNING THE MODEL
%% Obtain the divergence matrix
DE = cell(length(qN), length(qP), length(lambda), length(threshold));
DD = cell(length(qN), length(qP), length(lambda), length(threshold));
DN = cell(length(qN), length(qP), length(lambda), length(threshold));
iterations = 0;
totalCount = length(qN)*length(qP)*length(lambda)*length(threshold);
for i=1:length(qN)
    for j=1:length(qP)
        for k = 1:length(lambda)
            for l = 1:length(threshold)
                oneGridPtTime = tic;
                theta = [qN(i), qP(j), lambda(k), threshold(l)];
                [DE{i,j,k,l}, DD{i,j,k,l}, DN{i,j,k,l}] = divergence_endorsement_model(X,Z,Betas,leaders,TakeUp,EOmicronY,EOmicronM,DOmicronY,DOmicronM,NOmicronY,NOmicronM,Sec,theta,m,S,T,EmpRate,version,replicateOriginal);
                iterations = iterations+1;
                ['Done with ' num2str(iterations/totalCount*100) '% of the D(G,m) computations.']
                oneGridPtTimeelasped = toc(oneGridPtTime);
            end
        end
    end
    
end
toc;

% Process Data - construct data names
outputName = ['data_model_mom_' num2str(version) '']
save([outputName, ' ', timeVector, ' ','.mat'])

% 8. RUNNING THE AGGREGATOR
% Bootstrap?
bootstrap = 0;   % yes or no
if bootstrap==0,
    B = 1;
elseif bootstrap==1
    B = 1000;
end

%% Set up data
% Set up matrices
DivDD = zeros(length(qN), length(qP), length(lambda), length(threshold), G, m);
DivDE = zeros(length(qN), length(qP), length(lambda), length(threshold), G, m);
DivDN = zeros(length(qN), length(qP), length(lambda), length(threshold), G, m);

%% Put data into matrix

for i = 1:length(qN)
    for j = 1:length(qP)
        for k = 1:length(lambda)
            for l = 1:length(threshold)
                DDTotal(i, j, k, l, :, :) = DD{i, j, k, l};
                DETotal(i, j, k, l, :, :) = DE{i, j, k, l};
                DNTotal(i, j, k, l, :, :) = DN{i, j, k, l};
            end                
        end
    end
end

%Two Step Optimal Weights
twoStepOptimal = 0;

if twoStepOptimal == 1,
    theta = [qN_E, qP_E, lambda_E, threshold_E];
    [DE, DD, DN]= divergence_endorsement_model(X,Z,Betas,leaders,TakeUp,EOmicronY,EOmicronM,DOmicronY,DOmicronM,NOmicronY,NOmicronM,Sec,theta,m,S,T,EmpRate,version,replicateOriginal);
    AE = (DE'*DE)/43;
    WE = AE^(-1);
    
    theta = [qN_D, qP_D, lambda_D, threshold_D];
    [DE, DD, DN]= divergence_endorsement_model(X,Z,Betas,leaders,TakeUp,EOmicronY,EOmicronM,DOmicronY,DOmicronM,NOmicronY,NOmicronM,Sec,theta,m,S,T,EmpRate,version,replicateOriginal);
    AD = (DD'*DD)/43;
    WD = AD^(-1);
    
    theta = [qN_N, qP_N, lambda_N, threshold_N];
    [DE, DD, DN]= divergence_endorsement_model(X,Z,Betas,leaders,TakeUp,EOmicronY,EOmicronM,DOmicronY,DOmicronM,NOmicronY,NOmicronM,Sec,theta,m,S,T,EmpRate,version,replicateOriginal);
    AN = (DN'*DN)/43;
    WN = AN^(-1);
else
    WE = eye(m);
    WD = eye(m);
    WN = eye(m);
end

% Pre-allocation
QEndE = zeros(B,1);
QEndD = zeros(B,1);
QEndN = zeros(B,1);
TestEndEEndD = zeros(B,1);
TestEndDEndN = zeros(B,1);
TestEndEEndN = zeros(B,1);
importantparmsE = [];
importantparmsD = [];
importantparmsN = [];
valE = [];
valD = [];
valN = [];

% Aggregate
for b=1:B
    % Generate weights b
    if bootstrap==1,
        wt(b,:) = exprnd(1,G,1);
        wt(b,:) = wt(b,:)/mean(wt(b,:));
    elseif bootstrap==0
        wt(b,:) = 1/G*ones(1,G);
    end
    
    % For each model, generate the criterion function value for this Bootstrap to run
    
    % ENDORSEMENT MODEL
    for i=1:length(qN)
        for j=1:length(qP)
            for k=1:length(lambda)
                for l = 1:length(threshold)
                    % Compute the moment function
                    momFuncE(i,j,k,l,b,:) = (wt(b,:)*squeeze(DETotal(i,j,k,l,:,:)))'/G;
                    momFuncD(i,j,k,l,b,:) = (wt(b,:)*squeeze(DDTotal(i,j,k,l,:,:)))'/G;
                    momFuncN(i,j,k,l,b,:) = (wt(b,:)*squeeze(DNTotal(i,j,k,l,:,:)))'/G;              

                    % Criterion function
                    QEndorseE(i,j,k,l,b) = (squeeze(momFuncE(i,j,k,l,b,:)))'*WE*squeeze(momFuncE(i,j,k,l,b,:));
                    QEndorseD(i,j,k,l,b) = (squeeze(momFuncD(i,j,k,l,b,:)))'*WD*squeeze(momFuncD(i,j,k,l,b,:));
                    QEndorseN(i,j,k,l,b) = (squeeze(momFuncN(i,j,k,l,b,:)))'*WN*squeeze(momFuncN(i,j,k,l,b,:));
                end
            end
        end        
    end
    
    TempEndE = QEndorseE(:,:,:,:,b);
    TempEndD = QEndorseD(:,:,:,:,b);
    TempEndN = QEndorseN(:,:,:,:,b);
    
    [minAEndE,indEndE] = min(TempEndE(:));
    [x1EndE,x2EndE,x3EndE,x4EndE] = ind2sub(size(TempEndE),indEndE);
    [minAEndD,indEndD] = min(TempEndD(:));
    [x1EndD,x2EndD,x3EndD,x4EndD] = ind2sub(size(TempEndD),indEndD);
    [minAEndN,indEndN] = min(TempEndN(:));
    [x1EndN,x2EndN,x3EndN,x4EndN] = ind2sub(size(TempEndN),indEndN);
    
    importantparmsE = [importantparmsE; x1EndE, x2EndE, x3EndE, x4EndE];
    importantparmsD = [importantparmsD; x1EndD, x2EndD, x3EndD, x4EndD];
    importantparmsN = [importantparmsN; x1EndN, x2EndN, x3EndN, x4EndN];
    valE = [valE; qN(x1EndE), qP(x2EndE), lambda(x3EndE), threshold(x4EndE)];
    valD = [valD; qN(x1EndD), qP(x2EndD), lambda(x3EndD), threshold(x4EndD)];
    valN = [valN; qN(x1EndN), qP(x2EndN), lambda(x3EndN), threshold(x4EndN)];
    
    % Need to map back
    QEndE(b) = QEndorseE(x1EndE,x2EndE,x3EndE,x4EndE,b);
    QEndD(b) = QEndorseD(x1EndD,x2EndD,x3EndD,x4EndD,b);
    QEndN(b) = QEndorseN(x1EndN,x2EndN,x3EndN,x4EndN,b);
    
    % Tests
    TestEndEEndD(b) = sqrt(G)*(QEndE(b)-QEndD(b));
    TestEndDEndN(b) = sqrt(G)*(QEndD(b)-QEndN(b));
    TestEndEEndN(b) = sqrt(G)*(QEndE(b)-QEndN(b));
    
    ['Done with ' num2str(100*b/B) '% of the bootstraps']
    
end

outputNameNew = ['data_test_mom_' num2str(version) '']
save([outputNameNew, ' ', timeVector, ' ','.mat'])

%% The test
%{
CIendEEndD = [prctile(TestEndEEndD,5), prctile(TestEndEEndD,50), prctile(TestEndEEndD,95)]
CIendDEndN = [prctile(TestEndDEndN,5), prctile(TestEndDEndN,50), prctile(TestEndDEndN,95)]
CIendEEndN = [prctile(TestEndEEndN,5), prctile(TestEndEEndN,50), prctile(TestEndEEndN,95)]

[mean(valE); std(valE); median(valE)]
[mean(valD); std(valD); median(valD)]
[mean(valN); std(valN); median(valN)]

diffqNqPE = valE(:,1) - valE(:,2);
diffqNqPD = valD(:,1) - valD(:,2);
diffqNqPN = valN(:,1) - valN(:,2);

[mean(diffqNqPE); std(diffqNqPE)]
[mean(diffqNqPD); std(diffqNqPD)]
[mean(diffqNqPN); std(diffqNqPN)]
%}

