function FSMPrediction

% Evolutionary programming to optimize a finite state machine (FSM) to output a desired bit pattern 
% after seeing a specific input bit pattern.

% Individual format: [(output, next state), (output, next state), ...]
% We assume that the FSM begins in state 1.
% The first two integers refer to state 1 when the input is 0,
% the next two integers refer to state 1 when the input is 1,
% the next two integers refer to state 2 when the input is 0, ...
% So we use 4n integers to describe an FSM, where n is the number of states in the FSM, and
% the first state is number 1, and the last state is number n.

SequenceLength = 12;
RandStream.setDefaultStream(RandStream('mt19937ar', 'seed', sum(100*clock)));
TestInput = randi([0 1], 1, SequenceLength);
TestOutput = randi([0 1], 1, SequenceLength);
disp(['Test input =  ', num2str(TestInput)]);
disp(['Test output = ', num2str(TestInput)]);
    
beta = 1; gamma = 0; % EP parameters
NumGens = 100; % number of generations
PopSize = 5; % population size
NumStates = 4; % number of states in each individual of the evolving population
Pop = zeros(PopSize, 4*NumStates);
for i = 1 : PopSize
    for j = 1 : NumStates
        Pop(i, 4*(j-1)+1 : 4*j) = [randi([0 1]), randi([1 NumStates]), randi([0 1]), randi([1 NumStates])];
    end
end

MinCostArr = zeros(NumGens+1, 1);
for gen = 1 : NumGens
    % Evaluate each individual
    Cost = Evaluate(Pop, TestInput, TestOutput);
    MinCostArr(gen) = min(Cost);
    disp(['Min Cost = ', num2str(min(Cost))]);
    % Mutate
    CostNorm = 1 + (Cost - min(Cost)) / (max(Cost) - min(Cost)); % 1 <= CostNorm <= 2
    PopPrime = Pop + diag(sqrt(beta * CostNorm + gamma)) * randn(size(Pop));
    for i = 1 : PopSize
        for j = 1 : NumStates
            PopPrime(i, 4*(j-1)+1 : 4*j) = ...
                min(max(PopPrime(i, 4*(j-1)+1 : 4*j), [0 1 0 1]), [1 NumStates 1 NumStates]);
        end
    end
    PopPrime = round(PopPrime);
    % Evaluate the mutations
    CostPrime = Evaluate(PopPrime, TestInput, TestOutput);
    [~, ndx] = sort([Cost; CostPrime], 'ascend');
    PopNew = zeros(PopSize, 4*NumStates);
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
plot(0:NumGens, MinCostArr)
xlabel('Generation'), ylabel('Minimum Cost')
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Cost] = Evaluate(Pop, TestInput, TestOutput)
Cost = zeros(size(Pop, 1), 1);
for ind = 1 : size(Pop, 1)
    Cost(ind) = 0;
    % Evaluate the cost of each individual
    IndState = 1; % initial state of evolving individual
    for i = 1 : length(TestOutput)
        IndStateIndex = 4 * (IndState - 1) + 1;
        IndNextStateIndex = IndStateIndex + 2 * TestInput(i) + 1;
        IndState = Pop(ind, IndNextStateIndex);
        IndOutput = Pop(ind, IndNextStateIndex-1);
        Cost(ind) = Cost(ind) + abs(IndOutput - TestOutput(i));
    end
end
return