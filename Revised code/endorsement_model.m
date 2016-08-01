%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Return results of various stats, i.e., moments
% Four cases (with multiple moments) are programmed
%
% Originally programmed by Arun Chandrasekhar in Nov 2010
% Adapted by Hang Xiong in March 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% "parms" is renamed as "theta", the same as in other files.

function[infectedbeforeE,infectedbeforeD,infectedbeforeN,dynamicInfectionE,dynamicInfectionD,dynamicInfectionN,contagiousbeforeE,contagiousbeforeD,contagiousbeforeN] = endorsement_model(theta,Z,Betas,X,leaders,EOmicronY,EOmicronM,DOmicronY,DOmicronM,NOmicronY,NOmicronM,j,T,EmpRate,replicateOriginal)

%% Parameters: qN, qP, lambda, threshold (read them back from theta).
qN = theta(1);        % Probability non-taker transmits information.
qP = theta(2);        % Probability that a just-informed-taker transmits information.
lambda = theta(3);    % Coefficient in logit governing the probability of being a taker.
threshold = theta(4); % Threshold crossing which externalities take effects.

%% Pre-allocation
N = size(X,1);                  % Number of households in a village (graph)

% Eigenvector centrality
infectedE = false(N,1);         % Newly infected
infectedbeforeE = false(N,1);   % Infected previously
contagiousE = logical(leaders); % Newly informed/contagious
contagiousbeforeE = false(N,1); % Informed/contagious previously
transmissionHistE = false(N,N); % Who has transmitted to whom?
dynamicInfectionE = []; % A vector tracking the infection rate for the number of periods it takes place

% Degree centrality
infectedD = false(N,1);         % Newly infected
infectedbeforeD = false(N,1);   % Infected previously
contagiousD = logical(leaders); % Newly informed/contagious
contagiousbeforeD = false(N,1); % Informed/contagious previously
transmissionHistD = false(N,N); % Who has transmitted to whom?
dynamicInfectionD = []; % A vector tracking the infection rate for the number of periods it takes place

% Naive centrality
infectedN = false(N,1);         % Newly infected
infectedbeforeN = false(N,1);   % Infected previously
contagiousN = logical(leaders); % Newly informed/contagious
contagiousbeforeN = false(N,1); % Informed/contagious previously
transmissionHistN = false(N,N); % Who has transmitted to whom?
dynamicInfectionN = []; % A vector tracking the infection rate for the number of periods it takes place

%% Conduct the Simulation
x = rand(N,T);

% Eigenvector centrality
for t = 1:T
    % Step 1: Experience effect
    % construct regressor
    zE = infectedbeforeE * ones(1,N);
    regressorE = diag(EOmicronY*(transmissionHistE.*zE))./diag(EOmicronY*(transmissionHistE)); % fraction of adopting neighbors (F in the paper)
    regressorE(isnan(regressorE))=0;
    regressorE(isinf(regressorE))=0;
    LOGITprobE = 1./(1+exp(-([ones(N,1) Z]*Betas + regressorE*lambda))); % probability that an individual will adopt (p in the paper)

    infectedE = ((~contagiousbeforeE & contagiousE & x(:,t) < LOGITprobE) | infectedE);
    infectedbeforeE = (infectedbeforeE | infectedE);
    
    if(replicateOriginal==0)
        % Step 2: Externality effect
        infectedE = infectedE & (EOmicronM * infectedbeforeE < threshold); % Those who can borrow from their peers
        infectedbeforeE = (infectedbeforeE | infectedE);
    end
    
    contagiousbeforeE = (contagiousbeforeE | contagiousE);
    CE = sum(contagiousE);
    
    % Step 3: Information effect (get contagious)
    transmitPROBE = (contagiousE & infectedE)*qP + (contagiousE & ~infectedE)*qN; % Probability of transmission: Pr(.|infected)==qN, Pr(.|~infected)==qP. Individual (node) specific.
    contagionlikelihoodE = X.*(transmitPROBE*ones(1,N)); % A full NxN matrix of transmission probabilities
    t0E = rand(N,N);
    t1E = (contagionlikelihoodE > t0E); % transmission due to information effect 
    
    % zero out stuff because one can only transmit if contagious, and one
    % cannot be transmitted to unless they were not contagious before
    t1E(contagiousE,~contagiousbeforeE) =0; % This does not need to be done, but it makes no difference.
    t1E(~contagiousE,~contagiousbeforeE)=0; % If not contagious before, one cannot be transmitted.
    t1E(~contagiousE,contagiousbeforeE) =0; % Contagious before, not transmitted this time - no transmission
  
    t2E = t1E(contagiousE,~contagiousbeforeE); % which contagious folk transmit to previously uncontagious
    contagiousE(~contagiousbeforeE) = ((t2E'*ones(CE,1) > 0));
    
    transmissionHistE = (transmissionHistE | t1E);
    dynamicInfectionE = [dynamicInfectionE; sum(infectedbeforeE)/N];
end

% Degree centrality
for t = 1:T
    % Step 1: Experience effect
    % construct regressor
    zD = infectedbeforeD * ones(1,N);
    regressorD = diag(DOmicronY*(transmissionHistD.*zD))./diag(DOmicronY*(transmissionHistD)); % fraction of adopting neighbors (F in the paper)
    regressorD(isnan(regressorD))=0;
    regressorD(isinf(regressorD))=0;
    LOGITprobD = 1./(1+exp(-([ones(N,1) Z]*Betas + regressorD*lambda))); % probability that an individual will adopt (p in the paper)
 
    infectedD = (~contagiousbeforeD & contagiousD) & (x(:,t) < LOGITprobD | infectedD);
    infectedbeforeD = (infectedbeforeD | infectedD);
    
    if(replicateOriginal==0)
        % Step 2: Externality effect
        infectedD = infectedD & (DOmicronM * infectedbeforeD < threshold); % Those who can borrow from their peers
        infectedbeforeD = (infectedbeforeD | infectedD);
    end
    contagiousbeforeD = (contagiousbeforeD | contagiousD);
    CD = sum(contagiousD);
    
    % Step 3: Information effect (get contagious)
    transmitPROBE = (contagiousD & infectedD)*qP + (contagiousD & ~infectedD)*qN; % Probability of transmission: Pr(.|infected)==qN, Pr(.|~infected)==qP. Individual (node) specific.
    contagionlikelihoodD = X.*(transmitPROBE*ones(1,N)); % A full NxN matrix of transmission probabilities
    t0D = rand(N,N);
    t1D = (contagionlikelihoodD > t0D); % transmission due to information effect 
    
    % zero out stuff because one can only transmit if contagious, and one
    % cannot be transmitted to unless they were not contagious before
    t1D(contagiousD,~contagiousbeforeD) =0; % This does not need to be done, but it makes no difference.
    t1D(~contagiousD,~contagiousbeforeD)=0; % If not contagious before, one cannot be transmitted.
    t1D(~contagiousD,contagiousbeforeD) =0; % Contagious before, not transmitted this time - no transmission
  
    t2D = t1D(contagiousD,~contagiousbeforeD); % which contagious folk transmit to previously uncontagious
    contagiousD(~contagiousbeforeD) = ((t2D'*ones(CD,1) > 0));
    
    transmissionHistD = (transmissionHistD | t1D);
    dynamicInfectionD = [dynamicInfectionD; sum(infectedbeforeD)/N];
end

% Native centrality
for t = 1:T
    % Step 1: Experience effect
    % construct regressor
    zN = infectedbeforeN * ones(1,N);
    regressorN = diag(NOmicronY*(transmissionHistN.*zN))./diag(NOmicronY*(transmissionHistN)); % fraction of adopting neighbors (F in the paper)
    regressorN(isnan(regressorN))=0;
    regressorN(isinf(regressorN))=0;
    LOGITprobN = 1./(1+exp(-([ones(N,1) Z]*Betas + regressorN*lambda))); % probability that an individual will adopt (p in the paper)
 
    infectedN = (~contagiousbeforeN & contagiousN) & (x(:,t) < LOGITprobN | infectedN);
    infectedbeforeN = (infectedbeforeN | infectedN);
    
    if(replicateOriginal==0)
        % Step 2: Externality effect
        infectedN = infectedN & (NOmicronM * infectedbeforeN < threshold); % Those who can borrow from their peers
        infectedbeforeN = (infectedbeforeN | infectedN);
    end
    contagiousbeforeN = (contagiousbeforeN | contagiousN);
    CN = sum(contagiousN);
    
    % Step 3: Information effect (get contagious)
    transmitPROBE = (contagiousN & infectedN)*qP + (contagiousN & ~infectedN)*qN; % Probability of transmission: Pr(.|infected)==qN, Pr(.|~infected)==qP. Individual (node) specific.
    contagionlikelihoodN = X.*(transmitPROBE*ones(1,N)); % A full NxN matrix of transmission probabilities
    t0N = rand(N,N);
    t1N = (contagionlikelihoodN > t0N); % transmission due to information effect 
    
    % zero out stuff because one can only transmit if contagious, and one
    % cannot be transmitted to unless they were not contagious before
    t1N(contagiousN,~contagiousbeforeN) =0; % This does not need to be done, but it makes no difference.
    t1N(~contagiousN,~contagiousbeforeN)=0; % If not contagious before, one cannot be transmitted.
    t1N(~contagiousN,contagiousbeforeN) =0; % Contagious before, not transmitted this time - no transmission
  
    t2N = t1N(contagiousN,~contagiousbeforeN); % which contagious folk transmit to previously uncontagious
    contagiousN(~contagiousbeforeN) = ((t2N'*ones(CN,1) > 0));
    
    transmissionHistN = (transmissionHistN | t1N);
    dynamicInfectionN = [dynamicInfectionN; sum(infectedbeforeN)/N];
end

