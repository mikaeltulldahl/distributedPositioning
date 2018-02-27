function distributedPositioningSimulation()
global dt dim A Q
clc;
% rng(1) %set seed
dim = 2;
N = 20; %number of nodes
numberOfFixedNodes = 3;
numberOfParticles = 300;
timeSteps = 100;
noiseLevel = 0.01;
maxCost = 2.5*noiseLevel;
velocityDecay = 0.3;
listeningRange = 0.7;
plotRate = 2;
dt = 0.05;
A = [   1 0 dt 0;
    0 1 0 dt;
    0 0 (1-velocityDecay*dt) 0;
    0 0 0 (1-velocityDecay*dt)];
Q = diag(dt*[0.001 0.001 0.00003 0.00003]);


%init nodes
nodeStruct = struct('posTrue', {}, 'posEst', {}, 'posCov', {}...
    , 'velTrue', {}, 'velEst', {}, 'velCov', {}...
    , 'name', {}, 'timeStamp',{}, 'isFixed', {}...
    , 'knownNodes',{},'particles',{},'weights',{}...
    , 'dist',{}, 'distVar', {});
nodes = repmat(nodeStruct, N, 1 );
for n=1:N
    nodes(n).posTrue = rand(2,1);
    nodes(n).velTrue = zeros(2,1);
    nodes(n).name = n;
    nodes(n).timeStamp = 0;
    nodes(n).isFixed = false;
    nodes(n).knownNodes = nodeStruct;
    nodes(n).particles = [rand(2,numberOfParticles);0.05*rand(2,numberOfParticles)];
    nodes(n).weights = 1/numberOfParticles*ones(numberOfParticles, 1);
end

%set fixed nodes
fixed = randi(N,1,numberOfFixedNodes);
for n = 1:length(fixed)
    nodes(n).isFixed = fixed(n);
end

figure(2);
clf
grid on
semilogy(0, maxCost);
axis([0 timeSteps 0.1 1]);
hold on

for t = 1:timeSteps
    %update pos for each node
    for n = 1:N
        nodes(n).timeStamp = nodes(n).timeStamp + dt;
        if ~nodes(n).isFixed
            %update real pos
            newMean = A*[nodes(n).posTrue; nodes(n).velTrue];
            temp = mvnrnd(newMean',Q)';
            
            nodes(n).posTrue = temp(1:2);
            nodes(n).velTrue = temp(3:4);
        end
        nodes(n).knownNodes = []; %clear old knownNodes
        nodes(n) = propagate(nodes(n));
    end
    
    %each node broadcast position
    for n = 1:N
        %each node listen and records data + distance
        for m = 1:N
            if m ~= n && ~nodes(m).isFixed %can't hear itself, fixed nodes don't listen
                distTrue = norm(nodes(n).posTrue - nodes(m).posTrue);
                canHear = false;
                if distTrue < listeningRange
                    canHear = binornd(1,1-(distTrue/listeningRange)^2);
                end
                if canHear
                    noise = noiseLevel*randn;
                    distMeasured = distTrue + noise;
                    k = findNode(nodes(n), nodes(m).knownNodes);
                    newNode = nodes(n);
                    newNode.posTrue = [];
                    newNode.velTrue = [];
                    newNode.isFixed = [];
                    newNode.knownNodes = [];
                    newNode.particles = [];
                    newNode.weights = [];
                    newNode.dist = distMeasured;
                    newNode.distVar = 10*noiseLevel^2;
                    if k == 0 %insert
                        nodes(m).knownNodes = [nodes(m).knownNodes newNode];
                    else %overwrite
                        nodes(m).knownNodes(k) = newNode; %should not happen
                    end
                end % if canHear
            end% if m ~= n
        end % for m
    end %for n
    
    for n = 1:N
        if ~nodes(n).isFixed
            %update weights, normalize weights and resample using knownNodes
            nodes(n) = update(nodes(n));
        end
    end %for n
    
    clc
    errorSquareSum = 0;
    errorCounter = 0;
    for n = 1:N
        if ~nodes(n).isFixed
            error = norm(nodes(n).posEst - nodes(n).posTrue);
            errorSquareSum = errorSquareSum + error^2;
            errorCounter = errorCounter + 1;
        end
    end
    t
    normalizedErrorRMS = sqrt(errorSquareSum/errorCounter)/noiseLevel
%     for n = 1:N
%         if ~nodes(n).isFixed
%             cov = nodes(n).posCov
%             pos = nodes(n).posTrue
%             vel = nodes(n).velTrue
%         end
%     end
       
    figure(2);
    semilogy(t, normalizedErrorRMS,'.b');
    if mod((t+1),plotRate) == 0
        f = figure(1);
        clf
        axis([-0.3 1.3 -0.3 1.3]);
        grid on
        hold on
        for n = 1:N
            %plot where each node thinks it is and where it actually is
            figure(1);
            if nodes(n).isFixed
                plot(nodes(n).posTrue(1),nodes(n).posTrue(2),'om');
            else
                plot(nodes(n).posTrue(1),nodes(n).posTrue(2),'ob');
                plot(nodes(n).posEst(1),nodes(n).posEst(2),'xr');
%                 text(0.02+nodes(n).posEst(1),-0.02 + nodes(n).posEst(2),num2str(sqrt(nodes(n).posCov(1,1))*sqrt(nodes(n).posCov(2,2)),'%.2f'));
                text(0.02+nodes(n).posEst(1),nodes(n).posEst(2),num2str(nodes(n).name,'%u'));
                 plot_gaussian_ellipsoid(nodes(n).posEst, nodes(n).posCov, 2, 50, gca);
            end
        end% for n
    end %if plotThisRound
    drawnow;
    pause(0.1);
end
end

function node = propagate(node)
global dt A Q
if node.isFixed
    node.posEst = node.posTrue;
    node.posCov = 0.0001*eye(2);
    node.velEst = node.velTrue;
    node.velCov = zeros(2);
else
    for n = 1:size(node.particles,2)
        %propagate particles using motion model
        newParticleMean = A*node.particles(:,n);
        node.particles(:,n) = mvnrnd(newParticleMean',Q)';
    end
    %update estimate and cov using particles
    [node.posEst, node.posCov] = weightedMeanCov(node.particles(1:2,:), node.weights);
    [node.velEst, node.velCov] = weightedMeanCov(node.particles(3:4,:), node.weights);
end
node.timeStamp = node.timeStamp + dt;
end


function node = update(node)
if ~node.isFixed
    K = length(node.knownNodes);
    P = size(node.particles, 2);
    %update weights
    for k = 1:K
        for p = 1:P
            node.weights(p)= node.weights(p)*probabilityOfParticle(node.particles(1:2,p), node.knownNodes(k));
        end
    end
    
    %normalize weights
    node.weights = (1/sum(node.weights))*node.weights;
    
    %resample
    [node.particles, node.weights] = resampleTrashHistory(4, node.weights, node.particles);
end
end

function p = probabilityOfParticle(particlePos, node)
posDiff = particlePos - node.posEst;
direction = posDiff/norm(posDiff);
covariance = direction(1)*direction(2)*node.distVar + node.posCov(1,2);
varx = direction(1)^2*node.distVar + node.posCov(1,1);
vary = direction(2)^2*node.distVar + node.posCov(2,2);
positionCov = [ varx, covariance;
    covariance, vary];
p = mvnpdf(posDiff',(direction*node.dist)',positionCov);
end

function index = findNode(node, knownNodes)
K = length(knownNodes);
index = 0;
for k = 1:K
    if knownNodes(k).name == node.name
        index = k;
        break
    end
end
end

function [meanOut, covOut] = weightedMeanCov(samples, weights)
%https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_covariance
%assuming each row is a sampled random variable
numberOfVariables = size(samples,1);
numberOfSamples = size(samples,2);
weights = 1/sum(weights)*weights;
meanOut = samples*weights;
noMeanSamples = samples - repmat(meanOut, 1, numberOfSamples);
covOut = zeros(numberOfVariables,numberOfVariables);
for k = 1:numberOfSamples
    covOut = covOut + weights(k)*noMeanSamples(:,k)*noMeanSamples(:,k)';
end
end

function [ particlek, weightk ] = resampleTrashHistory(numberOfStates, weightk, particlek)
nrOfParticles = size(particlek, 2);
temp = sort(rand(nrOfParticles,1));
j=1;
weightSum=0;
resampledParticle = zeros(numberOfStates, nrOfParticles);
for i=1:nrOfParticles
    while weightSum <= temp(i)
        weightSum = weightSum + weightk(j);
        j=j+1;
    end
    resampledParticle(:,i) = particlek(:,j-1);
end
particlek = resampledParticle;
weightk = (1/nrOfParticles)*ones(nrOfParticles,1);
end
