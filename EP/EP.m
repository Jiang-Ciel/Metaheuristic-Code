function [MinCost] = EP(ProblemFunction, DisplayFlag, AdaptFlag, PopSize, InitPopSize, RandSeed, MaxGen)
% Evolutionary programming to optimize a general continuous function.
% INPUTS: ProblemFunction = the handle of the function that returns 
%                           the handles of the initialization, cost, and feasibility functions.
%         DisplayFlag = true or false, whether or not to display and plot results.
%         AdaptFlag = flag indicating whether or not to self-adapt mutation variance
%         PopSize = population size after initialization
%         InitPopSize = initial population size used to obtain initial population
%         RandSeed = random number seed
%         MaxGen = generation limit
% OUTPUT: MinCost = array of best solution, one element for each generation
if ~exist('ProblemFunction', 'var') || isempty(ProblemFunction)
    ProblemFunction = @Ackley;
end
if ~exist('DisplayFlag', 'var') || isempty(DisplayFlag)
    DisplayFlag = true;
end
if ~exist('AdaptFlag', 'var') || isempty(AdaptFlag)
    AdaptFlag = false;
end
if ~exist('PopSize', 'var') || isempty(PopSize)
    PopSize = 100;
end
if ~exist('InitPopSize', 'var') || isempty(InitPopSize)
    InitPopSize = PopSize;
end
InitPopSize = max(InitPopSize, PopSize);
if ~exist('RandSeed', 'var') || isempty(RandSeed)
    RandSeed = fix(sum(100*clock));
end
if ~exist('MaxGen', 'var') || isempty(MaxGen)
    MaxGen = 100;
end
% Initialization
OPTIONS.Maxgen = MaxGen;
OPTIONS.popsize = InitPopSize;
[OPTIONS, MinCost, AvgCost, Population, MinConstrViol, AvgConstrViol] = ...
    Init(DisplayFlag, ProblemFunction, OPTIONS, RandSeed);
Population = Population(1 : PopSize);
OPTIONS.popsize = PopSize;
OPTIONS.clearDups = false;
% EP parameters
beta = (max(OPTIONS.MaxDomain) - min(OPTIONS.MinDomain)) / 10;
gamma = 0; 
if AdaptFlag
    cParameter = 1;
    varMin = beta / 10;
    varMax = 10 * beta;
    for i = 1 : OPTIONS.popsize
        Population(i).var = varMin + (varMax - varMin) * rand;
    end
end
PopPrime = Population;
% Begin the optimization loop
for GenIndex = 1 : OPTIONS.Maxgen
    Cost = [Population.cost];
    CostNorm = 1 + (Cost - min(Cost)) / (max(Cost) - min(Cost)); % 1 <= CostNorm <= 2
    % EP Mutation
    for i = 1 : OPTIONS.popsize
        if AdaptFlag
            PopPrime(i).chrom = Population(i).chrom + sqrt(Population(i).var) * randn(1, OPTIONS.numVar);
            PopPrime(i).var = PopPrime(i).var + sqrt(cParameter * Population(i).var) * randn;
            PopPrime(i).var = max(PopPrime(i).var, varMin);
        else
            PopPrime(i).chrom = Population(i).chrom + sqrt(beta * CostNorm(i) + gamma) .* randn(1, OPTIONS.numVar);
        end
        % Make sure the mutated individuals are within the search domain
    	PopPrime(i).chrom = min( max( PopPrime(i).chrom, OPTIONS.MinDomain), OPTIONS.MaxDomain);
    end
    % Make sure the mutated population does not have duplicates. 
    if OPTIONS.clearDups
        PopPrime = ClearDups(PopPrime, OPTIONS);
    end
    % Calculate cost
    PopPrime = OPTIONS.CostFunction(PopPrime, OPTIONS);
    Population = PopSort([Population PopPrime], OPTIONS.ConstrMethod, GenIndex);
    Population = Population(1 : OPTIONS.popsize);
    [MinCost, AvgCost, MinConstrViol, AvgConstrViol] = ComputeCostAndConstrViol(Population, ...
        MinCost, AvgCost, MinConstrViol, AvgConstrViol, GenIndex, DisplayFlag);
end
Conclude(DisplayFlag, OPTIONS, Population, MinCost, AvgCost, MinConstrViol, AvgConstrViol);
return