function distributedPositioningSimulation()
clc;
clear;
rng(1354) %set seed
N = 20; %number of nodes
numberOfFixedNodes = 5;
timeSteps = 300;
noiseLevel = 0.02;
maxCost = 2.5*noiseLevel;
forgetRate = 1.1;
anchorAccuracy = 1*noiseLevel;
listeningRange = 0.8;
plotRate = 5;

%init nodes
knownNodes = struct('name',{},'posEst',{},'vel',{},'accuracy',{}, 'dist',{},'timeStamp',{});
nodes = repmat( struct('name', 0,'posEst',{},'posFit',{},'posTrue',{},'vel',{},'accuracy',{}, 'isFixed', false, 'knownNodes',knownNodes), N, 1 );
for n=1:N
    nodes(n).name = n;
    nodes(n).posEst = zeros(2,1);
    nodes(n).posFit = zeros(2,1);
    nodes(n).posTrue = rand(2,1);
    nodes(n).vel = zeros(2,1);
    nodes(n).accuracy = 1000;
    nodes(n).isFixed = false;
    if n <= numberOfFixedNodes
        nodes(n).isFixed = true;
        nodes(n).posEst = nodes(n).posTrue;
        nodes(n).accuracy = anchorAccuracy; %achor position accuracy
    end
    
end

figure(2);
clf
grid on
semilogy(0, maxCost);
axis([0 timeSteps 0.1 1]);
hold on

for t = 1:timeSteps
    plotThisRound = mod(t,plotRate) == 1;
    if plotThisRound
        figure(1);
        clf
        axis([-0.1 1.1 -0.1 1.1]);
        grid on
        hold on
    end
    
    %update pos for each node
    for n = (numberOfFixedNodes+1):N
        %update real pos
        
        %predict pos
        if ~nodes(n).isFixed
            timeDiff = 1;
            nodes(n).posEst = nodes(n).posEst + timeDiff*nodes(n).vel;
            nodes(n).vel = zeros(2,1);
            nodes(n).accuracy = forgetRate*nodes(n).accuracy;
        end
    end
    
    %each node broadcast position
    for n = 1:N
        %each node listen and records data + distance
        for m = 1:N
            if m ~=n %can't hear itself
                trueDistance = norm(nodes(n).posTrue - nodes(m).posTrue);
                canHear = false;
                if trueDistance < listeningRange
                    canHear = binornd(1,1-(trueDistance/listeningRange)^2);
                end
                if canHear
                    k = findNode(nodes(n).name, nodes(m).knownNodes);
                    overwrite = k == 0 || nodes(m).knownNodes(k).timeStamp < t;
                    if overwrite
                        newNode.name = nodes(n).name;
                        newNode.posEst = nodes(n).posEst;
                        newNode.vel = nodes(n).vel;
                        newNode.accuracy = nodes(n).accuracy;
                        newNode.dist = 0;
                        newNode.timeStamp = t;
                        if k == 0 %insert
                            nodes(m).knownNodes = [nodes(m).knownNodes newNode];
                            k = length(nodes(m).knownNodes);
                        else %overwrite
                            nodes(m).knownNodes(k) = newNode;
                        end
                    end
                    noise = noiseLevel*randn;
                    %perform distance measurement
                    nodes(m).knownNodes(k).dist = trueDistance + noise;
                end
            end
        end
        
        % update posEst, vel, accuracy of node n using knownNodes
        nodes(n) = update(nodes(n), maxCost);
        
        %plot where each node thinks it is and where it actually is
        if plotThisRound
            figure(1);
            if nodes(n).isFixed
                plot(nodes(n).posTrue(1),nodes(n).posTrue(2),'om');
            else
                plot(nodes(n).posTrue(1),nodes(n).posTrue(2),'ob');
                plot(nodes(n).posEst(1),nodes(n).posEst(2),'xr');
                plot(nodes(n).posFit(1),nodes(n).posFit(2),'xg');
                text(0.02+nodes(n).posEst(1),-0.02 + nodes(n).posEst(2),num2str(nodes(n).accuracy,'%.2f'));
                text(0.02+nodes(n).posEst(1),nodes(n).posEst(2),num2str(nodes(n).name,'%u'));
            end
        end
    end
    
    clc
    errorSquareSum = 0;
    for n = (numberOfFixedNodes+1):N
        error = norm(nodes(n).posEst - nodes(n).posTrue);
        errorSquareSum = errorSquareSum + error^2;
    end
    
    normalizedErrorRMS = sqrt(errorSquareSum/double(N - numberOfFixedNodes))/noiseLevel
    figure(2);
    semilogy(t, normalizedErrorRMS,'.b');
    drawnow;
end
end

function node = update(node, maxCost)
if ~node.isFixed
    % find best fit position for all distance measurments
    startEstimate = node.posEst;
    if node.accuracy >= maxCost %if old estimate is shit
        startEstimate = rand(2,1);
    end
    
    [node.posFit, cost] = fminsearch(@(pos)costFunction(pos, node.knownNodes),startEstimate);
    if cost < maxCost %if new estimate is OK
        if node.accuracy < maxCost %if old estimate was OK
            updateRatio = (atan(0.3*node.accuracy/cost*(pi/2))/(pi/2))^2;
            node.posEst = updateRatio*node.posFit + (1-updateRatio)*node.posEst;
            node.accuracy = updateRatio*cost + (1-updateRatio)*node.accuracy;
        else
            node.posEst = node.posFit;
            node.accuracy = cost;
        end
    end
end
end

function cost = costFunction(pos, knownNodes)
K = length(knownNodes);
sqrdError = 0;
weightSum = 0;
for k = 1:K
    measDist = knownNodes(k).dist;
    calcDist = norm(knownNodes(k).posEst - pos);
    error = calcDist - measDist;
    weight = (1/knownNodes(k).accuracy)^2;
    sqrdError = sqrdError + weight*error^2;
    weightSum = weightSum + weight;
end
cost = sqrt(sqrdError/weightSum);%rms
end

function index = findNode(name, knownNodes)
K = length(knownNodes);
index = 0;
for k = 1:K
    if knownNodes(k).name == name
        index = k;
        break
    end
end
end
