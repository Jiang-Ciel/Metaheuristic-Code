function [NumStates, MinCostArr] = PrimePrediction(GenLimit, SizePenalty, DisplayFlag)

% Evolutionary programming to optimize a finite state machine (FSM) to predict prime numbers.
% The input at each time step is the output of the previous time step.

% Individual format: [(start state), (output, next state), (output, next state), ...]
% The first integer is the starting state.
% The next two integers refer to state 1 when the input is 0,
% the next two integers refer to state 1 when the input is 1,
% the next two integers refer to state 2 when the input is 0, ...
% So we use 4n+1 integers to describe an FSM, where n is the number of states in the FSM.
% The first state is number 1, and the last state is number n.

if ~exist('GenLimit', 'var') || isempty(GenLimit)
    GenLimit = 100;
end
if ~exist('DisplayFlag', 'var') || isempty(DisplayFlag)
    DisplayFlag = true;
end
if ~exist('SizePenalty', 'var') || isempty(SizePenalty)
    SizePenalty = 0;
end

rng(round(sum(100*clock)));
SequenceLength = 100;
TestInput = isprime(1:SequenceLength-1);
TestOutput = isprime(2:SequenceLength);
    
PopSize = 20; % population size
NumStates = 4; % number of states in each individual of the initial population
Pop = zeros(PopSize, 4*NumStates+1);
for i = 1 : PopSize
    Pop(i, 1) = randi([1 NumStates]);
    for j = 1 : NumStates
        Pop(i, 4*(j-1)+2 : 4*j+1) = [randi([0 1]), randi([1 NumStates]), randi([0 1]), randi([1 NumStates])];
    end
end

MinCostArr = zeros(GenLimit+1, 1);
for gen = 1 : GenLimit
    % Evaluate each individual
    [Cost, ModCost] = Evaluate(Pop, TestInput, TestOutput, SizePenalty);
    MinCostArr(gen) = min(Cost);
    if DisplayFlag
        disp(['Min Cost = ', num2str(min(Cost))]); 
    end
    % Mutate
    PopPrime = Pop;
    for i = 1 : PopSize
        NumStates = CalcNumStates(PopPrime(i,:));
        MutationType = randi([1 5]);
        if MutationType == 1
            % add a state
            NewState = [randi([0 1]), randi([1 NumStates]), randi([0 1]), randi([1 NumStates])]; 
            NumColsBeforeAdding = size(PopPrime, 2);
            PopPrime(i, 4*NumStates+2 : 4*(NumStates+1)+1) = NewState;
            if size(PopPrime, 2) > NumColsBeforeAdding
                PopPrime(1:i-1, end-3:end) = inf;
                PopPrime(i+1:end, end-3:end) = inf;
            end
        elseif MutationType == 2
            % delete a state if there is more than one state
            if NumStates > 1
                DeleteState = randi([1 NumStates]);
                PopPrime(i, 4*(DeleteState-1)+2 : end-4) = PopPrime(i, 4*DeleteState+2 : end);
                PopPrime(i, 4*NumStates-2 : 4*NumStates+1) = inf;
                if PopPrime(i, 1) == DeleteState
                    PopPrime(i, 1) = randi([1, NumStates-1]);
                elseif PopPrime(i, 1) > DeleteState
                    PopPrime(i, 1) = PopPrime(i, 1) - 1;
                end
                for j = 1 : NumStates-1
                    if PopPrime(i, 4*j-1) == DeleteState
                        PopPrime(i, 4*j-1) = randi([1, NumStates-1]);
                    elseif PopPrime(i, 4*j-1) > DeleteState
                        PopPrime(i, 4*j-1) = PopPrime(i, 4*j-1) - 1;
                    end
                    if PopPrime(i, 4*j+1) == DeleteState
                        PopPrime(i, 4*j+1) = randi([1, NumStates-1]);
                    elseif PopPrime(i, 4*j+1) > DeleteState
                        PopPrime(i, 4*j+1) = PopPrime(i, 4*j+1) - 1;
                    end
                end
            end
        elseif MutationType == 3
            % change an output
            PopPrime(i, 4*randi([0 NumStates-1])+2+2*randi([0 1])) = randi([0 1]);
        elseif MutationType == 4
            % change a state transition
            PopPrime(i, 4*randi([0 NumStates-1])+1+2*randi([1 2])) = randi([1 NumStates]);
        else
            % change the starting state
            PopPrime(i, 1) = randi([1 NumStates]);
        end
    end
    % Evaluate the mutations
    [CostPrime, ModCostPrime] = Evaluate(PopPrime, TestInput, TestOutput, SizePenalty);
    [~, ndx] = sort([ModCost; ModCostPrime], 'ascend');
    PopNew = inf(PopSize, max(size(Pop,2), size(PopPrime,2)));
    for i = 1 : PopSize
        if ndx(i) <= PopSize
            PopNew(i, 1:size(Pop,2)) = Pop(ndx(i), :);
        else
            PopNew(i, 1:size(PopPrime,2)) = PopPrime(ndx(i)-PopSize, :);
        end
    end
    Pop = PopNew;
end
MinCostArr(gen+1) = min([Cost; CostPrime]);
NumStates = CalcNumStates(Pop(1,:));
if DisplayFlag
    disp(['Best FSM has ', num2str(NumStates), ' states: ', num2str(Pop(1, :))]);
    close all; SetPlotOptions
    plot(0:GenLimit, MinCostArr)
    xlabel('Generation'), ylabel('Minimum Cost')
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Cost, ModCost] = Evaluate(Pop, TestInput, TestOutput, SizePenalty)
Cost = zeros(size(Pop, 1), 1);
ModCost = zeros(size(Pop, 1), 1);
for ind = 1 : size(Pop, 1)
    Cost(ind) = 0;
    % Evaluate the cost of each individual
    IndState = Pop(ind, 1); % initial state of evolving individual
    for i = 1 : length(TestOutput)
        IndStateIndex = 4 * (IndState - 1) + 2;
        IndNextStateIndex = IndStateIndex + 2 * TestInput(i) + 1;
        IndState = Pop(ind, IndNextStateIndex);
        IndOutput = Pop(ind, IndNextStateIndex-1);
        Cost(ind) = Cost(ind) + abs(IndOutput - TestOutput(i));
    end
    ModCost(ind) = SizePenalty * CalcNumStates(Pop(ind,:)) + Cost(ind);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NumStates] = CalcNumStates(FSMVector)
NumStates = find(FSMVector==inf, 1);
if isempty(NumStates)
    NumStates = (length(FSMVector) - 1) / 4;
else
    NumStates = (NumStates - 2) / 4;
end
return