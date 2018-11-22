function [FoodEaten, Moves, BestFSM] = SanteFe32(NumStates, DisplayFlag)

% Evolutionary programming to optimize a finite state machine (FSM) for navigating an
% artificial ant on the 32 x 32 Sante Fe trail. See Section 3.3.2 in "Genetic Programming," 
% by John Koza, for the problem statement.

% Individual format: [(output, next state), (output, next state), ...]
% The first two integers refer to state 1 when the input is not 1 (no food sensed),
% the next two integers refer to state 1 when the input is 1 (food is sensed),
% the next two integers refer to state 2 when the input is not 1, ...
% So there are 4n integers to describe an FSM, where n is the number of states in the FSM.
% The first state is state number 1, and the last state is state number n.
% Outputs are coded as 0 for move forward, 1 for turn right, and 2 for turn left.
% If an element in a strategy is equal to infinity, then there are no more valid states in that strategy.

MaxMoves = 500;
FirstPass = true;
if ~exist('NumStates', 'var') || isempty(NumStates)
    NumStates = 5;
end
if ~exist('DisplayFlag', 'var') || isempty(DisplayFlag)
    DisplayFlag = true;
end
beta = 1; gamma = 0; % EP parameters
NumGens = 100; % number of generations
PopSize = 100; % population size
% Initialize the population of FSMs
Pop = zeros(PopSize, 4*NumStates);
PopNew = Pop;
for i = 1 : PopSize
    for j = 1 : NumStates
        Pop(i, 4*j-3 : 4*j) = [randi([0 2]), randi([1 NumStates]), randi([0 2]), randi([1 NumStates])];
    end
end
rng(round(sum(100*clock))); % initialize the random number generator
MinCostArr = zeros(1, NumGens+1);
AveCostArr = zeros(1, NumGens+1);
for gen = 0 : NumGens
    % Evaluate each individual
    [FoodEaten, Moves, FirstPass] = Evaluate(Pop, MaxMoves, FirstPass, DisplayFlag);
    Cost = Moves - min(Moves) + max(FoodEaten) - FoodEaten;
    MinCostArr(gen+1) = max(FoodEaten);
    AveCostArr(gen+1) = mean(FoodEaten);
    if DisplayFlag
        disp(['Generation ', num2str(gen), ', most food eaten = ', num2str(max(FoodEaten))])
    end
    if (gen == NumGens), break, end
    % Mutate
    CostNorm = 1 + (Cost - min(Cost)) / (max(Cost) - min(Cost)); % 1 <= CostNorm <= 2
    PopPrime = Pop + diag(sqrt(beta * CostNorm + gamma)) * randn(size(Pop));
    for i = 1 : PopSize
        for j = 1 : NumStates
            PopPrime(i, 4*j-3 : 4*j) = ...
                min(max(PopPrime(i, 4*j-3 : 4*j), [0 1 0 1]), [2 NumStates 2 NumStates]);
        end
    end
    PopPrime = round(PopPrime);
    % Evaluate the mutations
    [FoodEatenPrime, MovesPrime, FirstPass] = Evaluate(PopPrime, MaxMoves, FirstPass, DisplayFlag);
    CostPrime = MovesPrime - min([Moves; MovesPrime]) + max([FoodEaten; FoodEatenPrime]) - FoodEatenPrime;
    Cost = Moves - min([Moves; MovesPrime]) + max([FoodEaten; FoodEatenPrime]) - FoodEaten;
    [~, ndx] = sort([Cost; CostPrime], 'ascend');
    for i = 1 : PopSize
        if ndx(i) <= PopSize
            % Add one of the parent individuals to the new population 
            PopNew(i, :) = Pop(ndx(i), :);
        else
            % Add one of the child individuals to the new population
            PopNew(i, :) = PopPrime(ndx(i)-PopSize, :);
        end
    end
    Pop = PopNew;
end
if DisplayFlag
    disp(['Best FSM = ', num2str(Pop(1, :))]);
    figure, plot(0:NumGens,MinCostArr,'r-', 0:NumGens,AveCostArr,'b--')
    xlabel('Generation'), ylabel('Fitness')
    legend('Maximum Fitness', 'Average Fitness')
end
if ndx(1) <= PopSize
    FoodEaten = FoodEaten(ndx(1));
    Moves = Moves(ndx(1));
else
    FoodEaten = FoodEatenPrime(ndx(1)-PopSize);
    Moves = MovesPrime(ndx(1)-PopSize);
end
BestFSM = Pop(1, :);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FoodEaten, Moves, FirstPass] = Evaluate(Pop, MaxMoves, FirstPass, DisplayFlag)
FoodEaten = zeros(size(Pop, 1), 1);
Moves = zeros(size(Pop, 1), 1);
% Evaluate the performance of each finite state machine represented in the Pop array
for ind = 1 : size(Pop, 1)
    [Food, FirstPass] = CreateFoodArray(FirstPass, DisplayFlag); % Initialize the food array
    Rows = size(Food, 1); % number of rows in toroidal grid
    Cols = size(Food, 2); % number of columns in toroidal grid
    FoodAvailable = sum(Food(:) == 1); % how much food is in the grid
    MoveNum = 0; % how many moves the ant has taken
    FoodEaten(ind) = 0; % how much food has been eaten so far
    State = 1; % initial state
    Square = [1, 1]; % initial ant location [row, column]
    Direction = 0; % 0 means facing right, 1 means facing up, 2 means facing left, 3 means facing down
    while (MoveNum < MaxMoves) && (FoodEaten(ind) < FoodAvailable)
        % Calculate NextSquare, the [row column] indices of the cell which the ant is facing
        if Direction == 0
            NextSquare = Square + [0 1];
        elseif Direction == 1
            NextSquare = Square + [1 0];
        elseif Direction == 2
            NextSquare = Square - [0 1];
        else
            NextSquare = Square - [1 0];
        end
        if NextSquare(1) > Rows, NextSquare(1) = 1; end
        if NextSquare(1) < 1, NextSquare(1) = Rows; end
        if NextSquare(2) > Cols, NextSquare(2) = 1; end
        if NextSquare(2) < 1, NextSquare(2) = Cols; end
        Ndx = 4 * State - 3; % index into the first element of Pop(ind) that corresponds to the current state
        FoodIndicator = Food(NextSquare(1), NextSquare(2)); % does the ant sense food?
        if FoodIndicator == 1
            % Food ahead
            Action = Pop(ind, Ndx+2);
            State = Pop(ind, Ndx+3);
        else
            % No food ahead
            Action = Pop(ind, Ndx);
            State = Pop(ind, Ndx+1);
        end
        if Action == 0
            % move forward
            Square = NextSquare;
            if FoodIndicator == 1
                FoodEaten(ind) = FoodEaten(ind) + 1;
                Food(Square(1), Square(2)) = 0;
            end
        elseif Action == 1
            % turn right
            Direction = Direction - 1;
            if Direction < 0, Direction = 3; end
        else
            % turn left
            Direction = Direction + 1;
            if Direction > 3, Direction = 0; end
        end
        MoveNum = MoveNum + 1;
    end
    Moves(ind) = MoveNum;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Food, FirstPass] = CreateFoodArray(FirstPass, DisplayFlag)
Food = zeros(32, 32); % 0 = no food, 1 = food, 2 = no food but optimal path
Food(1, 2:4) = 1;
Food(2:6, 4) = 1;
Food(6, 5:7) = 1;
Food(6, 8) = 2;
Food(6, 9:13) = 1;
Food(7:11, 13) = 1;
Food(12, 13) = 2;
Food(13:16, 13) = 1;
Food(17:18, 13) = 2;
Food(19:25, 13) = 1;
Food(26, 13) = 2;
Food(26, 8:12) = 1;
Food(26, 6:7) = 2;
Food(26, 4:5) = 1;
Food(26, 2:3) = 2;
Food(27:30, 2) = 1;
Food(31:32, 2) = 2;
Food(32, 3:6) = 1;
Food(32, 7:9) = 2;
Food(30:31, 9) = 1;
Food(29, 9) = 2;
Food(29, 10:16) = 1;
Food(29, 17:18) = 2;
Food(26:28, 18) = 1;
Food(24:25, 18) = 2;
Food(20:23, 18) = 1;
Food(17:19, 18) = 2;
Food(17, 19) = 1;
Food(17, 20:22) = 2;
Food(15:16, 22) = 1;
Food(13:14, 22) = 2;
Food(9:12, 22) = 1;
Food(7:8, 22) = 2;
Food(7, 23:24) = 1;
Food(7, 25:26) = 2;
Food(5:6, 26) = 1;
Food(4, 26) = 2;
Food(4, 27:29) = 1;
Food(4, 30:31) = 2;
Food(5:6, 31) = 1;
Food(7, 31) = 2;
Food(8, 31) = 1;
Food(9:10, 31) = 2;
Food(11, 31) = 1;
Food(12:13, 31) = 2;
Food(14, 31) = 1;
Food(15:16, 31) = 2;
Food(16, 28:30) = 1;
Food(16, 25:27) = 2;
Food(17, 25) = 1;
Food(18:20, 25) = 2;
Food(20, 26) = 1;
Food(20, 27:29) = 2;
Food(21, 29) = 1;
Food(22:24, 29) = 2;
Food(24, 28) = 1;
Food(24, 25:27) = 2;
Food(25, 25) = 1;
if FirstPass
    FirstPass = false;
    if DisplayFlag
        SetPlotOptions
        colormap([1 1 1; 0 0 0; .5 .5 .5])
        FoodPlot = [Food, zeros(32, 1)];
        FoodPlot = [FoodPlot; zeros(1, 33)];
        pcolor(FoodPlot)
    end
end
return
