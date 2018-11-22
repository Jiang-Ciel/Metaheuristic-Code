function Prisoner(SelfEvolve)

% Evolutionary programming to optimize a finite state machine (FSM) for the prisoner's dilemma problem.
% INPUT: SelfEvolve = true if the population evolves by playing against each other,
%                     false if the population evolves by playing against fixed opponents.

% Individual format: [initial output, [output, next state], [output, next state], ...]
% The first integer refers to the initial output. We assume that the FSM transitions to state 1 after
% its initial output.
% The next two integers refer to state 1 when the input is cooperate,
% the next two integers refer to state 1 when the input is defect,
% the next two integers refer to state 2 when the input is cooperate, ...
% So there are 4n+1 integers to describe an FSM, where n is the number of states in the FSM, where
% the first state is number 1, and the last state is number n.
% Outputs are coded as 0 for cooperate, and 1 for defect.
% If an element in a strategy is infinity, then there are no more valid states in that strategy.
% Assume that all strategies begin in state #1.
% Test strategies:
%  Strategy = [0, 0 1 0 1, inf*ones(1,12); % always cooperate
%              1, 1 1 1 1, inf*ones(1,12); % always defect
%              0, 0 1 1 1, inf*ones(1,12); % tit-for-tat
%              0, 0 1 0 2, 0 1 1 2, inf*ones(1,8); % tit-for-two-tats
%              0, 0 1 1 2, 1 2 1 2, inf*ones(1,8); % grim
%              0, 0 1 1 2, 1 3 1 2, 1 4 1 2, 0 1 1 2]; % punish

if ~exist('SelfEvolve', 'var') || isempty(SelfEvolve)
    SelfEvolve = false;
end

if ~SelfEvolve
    NumStates = 4;
    NumStrategies = 4;
    Strategy = zeros(NumStrategies, 4*NumStates+1);
    Strategy(NumStrategies, 1) = randi([0 1]);
    for i = 1 : NumStrategies
        for j = 1 : NumStates
            Strategy(i, 4*(j-1)+2 : 4*j+1) = [randi([0 1]), randi([1 NumStates]), randi([0 1]), randi([1 NumStates])];
        end
        disp(['Opponent FSM(', num2str(i), ') = ', num2str(Strategy(i,:))]);
    end
end

beta = 1; gamma = 0; % EP parameters
NumGames = 10; % number of games
NumGens = 100; % number of generations
PopSize = 5; % population size
NumStates = 2; % number of states in each individual of the evolving population
Pop = zeros(PopSize, 4*NumStates+1);
for i = 1 : PopSize
    Pop(i, 1) = randi([0 1]);
    for j = 1 : NumStates
        Pop(i, 4*(j-1)+2 : 4*j+1) = [randi([0 1]), randi([1 NumStates]), randi([0 1]), randi([1 NumStates])];
    end
end
CostMatrix = [1 10; 0 5];

rng(round(sum(100*clock)));
MinCostArr = zeros(1, NumGens+1);
AveCostArr = zeros(1, NumGens+1);
for gen = 1 : NumGens
    % Evaluate each individual
    if SelfEvolve
        Strategy = Pop;
    end
    Cost = Evaluate(Pop, Strategy, CostMatrix, NumGames);
    MinCostArr(gen) = min(Cost);
    AveCostArr(gen) = mean(Cost);
    % Mutate
    CostNorm = 1 + (Cost - min(Cost)) / (max(Cost) - min(Cost)); % 1 <= CostNorm <= 2
    PopPrime = Pop + diag(sqrt(beta * CostNorm + gamma)) * randn(size(Pop));
    PopPrime(:, 1) = min(max(PopPrime(:, 1), 0), 1);
    for i = 1 : PopSize
        for j = 1 : NumStates
            PopPrime(i, 4*(j-1)+2 : 4*j+1) = ...
                min(max(PopPrime(i, 4*(j-1)+2 : 4*j+1), [0 1 0 1]), [1 NumStates 1 NumStates]);
        end
    end
    PopPrime = round(PopPrime);
    % Evaluate the mutations
    CostPrime = Evaluate(PopPrime, Strategy, CostMatrix, NumGames);
    [SortedCost, ndx] = sort([Cost; CostPrime], 'ascend');
    PopNew = zeros(PopSize, 4*NumStates+1);
    for i = 1 : PopSize
        if ndx(i) <= PopSize
            PopNew(i, :) = Pop(ndx(i), :);
        else
            PopNew(i, :) = PopPrime(ndx(i)-PopSize, :);
        end
    end
    Pop = PopNew;
end
disp(['Best FSM = ', num2str(Pop(1, :))]);
close all; SetPlotOptions
MinCostArr(gen+1) = min([Cost; CostPrime]);
AveCostArr(gen+1) = mean(SortedCost(ndx(1:PopSize)));
figure, plot(0:NumGens, MinCostArr)
xlabel('Generation'), ylabel('Minimum Cost')
figure, plot(0:NumGens, AveCostArr)
xlabel('Generation'), ylabel('Average Cost')
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Cost] = Evaluate(Pop, Strategy, CostMatrix, NumGames)
Cost = zeros(size(Pop, 1), 1);
for ind = 1 : size(Pop, 1)
    Cost(ind) = 0;
    % Evaluate the total cost of each individual after playing all opponent strategies
    for s = 1 : size(Strategy, 1)
        IndState = 1; % initial state of evolving individual
        IndMove = Pop(ind, 1); % initial move of evolving individual
        OppState = 1; % initial state of opponent
        OppMove = Strategy(s, 1); % initial move of opponent
        Cost(ind) = Cost(ind) + CostMatrix(IndMove+1, OppMove+1);
        for i = 1 : NumGames-1
            IndStateIndex = 4 * (IndState - 1) + 2;
            IndNextStateIndex = IndStateIndex + 2 * OppMove + 1;
            IndState = Pop(ind, IndNextStateIndex);
            OppStateIndex = 4 * (OppState - 1) + 2;
            OppNextStateIndex = OppStateIndex + 2 * IndMove + 1;
            OppState = Strategy(s, OppNextStateIndex);
            IndMove = Pop(ind, IndNextStateIndex-1);
            OppMove = Strategy(s, OppNextStateIndex-1);
            Cost(ind) = Cost(ind) + CostMatrix(IndMove+1, OppMove+1);
        end
    end
end
return