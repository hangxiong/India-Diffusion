%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Return the differences between empirical and simulated moments
% Call "endorsement_model.m" for running the simulations
% and "moments.m" for calculating the results of moments
%
% Originally programmed by Arun Chandrasekhar in May 2011
% Adapted by Hang Xiong in March 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [DE, DD, DN] = divergence_endorsement_model(X,Z,Betas,leaders,TakeUp,EOmicronY,EOmicronM,DOmicronY,DOmicronM,NOmicronY,NOmicronM,Sec,theta,m,S,T,EmpRate,version,replicateOriginal)


%% Parameters
G = length(X); % Number of graphs/villages in question.

%% Computation of the vector of divergences across all the moments
EmpiricalMoments = zeros(G,m);
MeanSimulatedMomentsE = zeros(G,m);
MeanSimulatedMomentsD = zeros(G,m);
MeanSimulatedMomentsN = zeros(G,m);
DE = zeros(G,m);
DD = zeros(G,m);
DN = zeros(G,m);

for g=1:G
    % Compute moments - G x m object
    EmpiricalMoments(g,:) = moments(X{g},leaders{g},TakeUp{g},Sec{g},g,version);
    % Compute simulated moments
    SimulatedMomentsE = zeros(S,m);
    SimulatedMomentsD = zeros(S,m);
    SimulatedMomentsN = zeros(S,m);
    for s=1:S
        oneSimTime = tic;
        [InfectedSimsE, InfectedSimsD, InfectedSimsN] = endorsement_model(theta,Z{g},Betas,X{g},leaders{g},EOmicronY{g},EOmicronM{g},DOmicronY{g},DOmicronM{g},NOmicronY{g},NOmicronM{g},g,T(g),EmpRate(g),replicateOriginal);
        SimulatedMomentsE(s,:) = moments(X{g},leaders{g},InfectedSimsE,Sec{g},g,version);
        SimulatedMomentsD(s,:) = moments(X{g},leaders{g},InfectedSimsD,Sec{g},g,version);
        SimulatedMomentsN(s,:) = moments(X{g},leaders{g},InfectedSimsN,Sec{g},g,version);
        toc(oneSimTime)
        [s g theta]
        ['Done with ' num2str(g/G*100) '% of the graphs for THIS parameter value']
    end

    % Compute the mean simulated moment - a G x m object
    MeanSimulatedMomentsE(g,:) = mean(SimulatedMomentsE,1);
    MeanSimulatedMomentsD(g,:) = mean(SimulatedMomentsD,1);
    MeanSimulatedMomentsN(g,:) = mean(SimulatedMomentsN,1);
    
    DE(g,:) = MeanSimulatedMomentsE(g,:) - EmpiricalMoments(g,:);
    DD(g,:) = MeanSimulatedMomentsD(g,:) - EmpiricalMoments(g,:);
    DN(g,:) = MeanSimulatedMomentsN(g,:) - EmpiricalMoments(g,:);
end

