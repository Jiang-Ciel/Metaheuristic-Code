function [MinCost] = BBO(ProblemFunction, DisplayFlag, LinearFlag, Blending, ConstrMethod, RandSeed, popsize, Maxgen)
% Biogeography-based optimization (BBO) software for minimizing a continuous function

% INPUTS: ProblemFunction = the handle of the function that returns 
%                           the handles of the initialization and cost functions.
%         DisplayFlag = true or false, whether or not to display and plot results.
%         LinearFlag = true to use linear migration models, false to use sinusoidal migration models
%         Blending = migration blend parameter, between 0 and 1; 0 means standard BBO
%         ConstrMethod = constraint-handling method for constrained optimization problems -
%                        see PopSort.m for description and algorithms of constraint-handling methods
%         RandSeed = random number seed
%         popsize = population size
%         Maxgen = generation limit
% OUTPUT: MinCost = array of best solution, one element for each generation
 
if ~exist('ProblemFunction', 'var') || isempty(ProblemFunction)
    ProblemFunction = @Ackley;
end
if ~exist('DisplayFlag', 'var') || isempty(DisplayFlag)
    DisplayFlag = true;
end
if ~exist('LinearFlag', 'var') || isempty(LinearFlag)
    LinearFlag = true;
end
if ~exist('Blending', 'var') || isempty(Blending)
    Blending = 0;
end
if exist('ConstrMethod', 'var') && ~isempty(ConstrMethod)
    OPTIONS.ConstrMethod = ConstrMethod;
else
    OPTIONS = [];
end
if ~exist('RandSeed', 'var') || isempty(RandSeed)
    RandSeed = [];
end
if ~exist('popsize', 'var') || isempty(popsize)
    popsize = 100;
end
OPTIONS.popsize = popsize;
if ~exist('Maxgen', 'var') || isempty(Maxgen)
    Maxgen = 100;
end
OPTIONS.Maxgen = Maxgen;

% Initialization
[OPTIONS, MinCost, AvgCost, Population, MinConstrViol, AvgConstrViol] = ...
    Init(DisplayFlag, ProblemFunction, OPTIONS, RandSeed);

% Compute immigration rate and emigration rate for each individual.
% lambda(i) is the immigration rate for habitat i.
% mu(i) is the emigration rate for habitat i.
[lambda, mu] = GetLambdaMu(length(Population), LinearFlag);

% Begin the optimization loop
for GenIndex = 1 : OPTIONS.Maxgen
    % Save the best habitats in a temporary array.
    % For unconstrained problems, save the best individuals.
    % For constrained problems, save the best feasible individuals, followed by the best infeasible individuals.
    ElitePop = Population(1 : OPTIONS.Keep);
    if isfield(Population, 'G')
        Feasibles = Population([Population.G] == 0);
        numFeas = length(Feasibles);
        for i = 1 : min(OPTIONS.Keep, numFeas);
            ElitePop(i) = Feasibles(i);
        end
        if OPTIONS.Keep > numFeas
            Infeasibles = Population([Population.G] ~= 0);
            for i = numFeas+1 : OPTIONS.Keep
                ElitePop(i) = Infeasibles(i-numFeas);
            end
        end
    end
    % Begin the migration loop
    TempPop = Population;
    for k = 1 : length(Population)
        % Probabilistically input new information into habitat i
        for j = 1 : OPTIONS.numVar
            if rand < lambda(k)
                % Pick a habitat from which to obtain a feature
                RandomNum = rand * sum(mu);
                Select = mu(1);
                SelectIndex = 1;
                while (RandomNum > Select) && (SelectIndex < OPTIONS.popsize)
                    SelectIndex = SelectIndex + 1;
                    Select = Select + mu(SelectIndex);
                end
                TempPop(k).chrom(j) = ...
                    Blending * Population(k).chrom(j) + (1 - Blending) * Population(SelectIndex).chrom(j);
            end
        end
    end
    % Mutation
    for k = 1 : length(Population)
        for j = 1 : OPTIONS.numVar
            if OPTIONS.pmutate > rand
                TempPop(k).chrom(j) = OPTIONS.MinDomain(j) + (OPTIONS.MaxDomain(j) - OPTIONS.MinDomain(j)) * rand;
            end
        end
    end
    % Replace the habitats with their new versions
    Population = TempPop;
    % Make sure the population does not have duplicates
    if OPTIONS.clearDups
        Population = ClearDups(Population, OPTIONS);
    end
    % Calculate cost
    Population = OPTIONS.CostFunction(Population, OPTIONS);
    % Replace the worst with the previous generation's elites
    Population = PopSort([Population, ElitePop], OPTIONS.ConstrMethod, GenIndex);
    Population = Population(1 : OPTIONS.popsize);
    [MinCost, AvgCost, MinConstrViol, AvgConstrViol] = ComputeCostAndConstrViol(Population, ...
        MinCost, AvgCost, MinConstrViol, AvgConstrViol, GenIndex, DisplayFlag);
end
Conclude(DisplayFlag, OPTIONS, Population, MinCost, AvgCost, MinConstrViol, AvgConstrViol);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lambda, mu] = GetLambdaMu(N, Linear)
% Compute immigration rate and extinction rate for each species count.
% lambda(i) is the immigration rate for individual i.
% mu(i) is the extinction rate for individual i.
% This routine assumes the population is sorted from most fit to least fit.
iArray = 1 : N;
if Linear
    mu = (N - iArray) / N; % linear migration curves
else
    mu = (1 + cos(iArray*pi/N)) / 2; % sinusoidal migration curves
end
lambda = 1 - mu;
return