function EPMonte(ProblemFunction)
% Monte Carlo simulation of EP to compare adaptation vs. no adaptation
nMonte = 20;
if ~exist('ProblemFunction', 'var') || isempty(ProblemFunction)
    ProblemFunction = @Ackley;
end
DisplayFlag = false;
PopSize = 100;
InitPopSize = 100;
MaxGen = 50;
MinCostNoAdapt = zeros(MaxGen+1, nMonte);
MinCostAdapt = zeros(MaxGen+1, nMonte);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    RandSeed = round(sum(clock));
    MinCostNoAdapt(:, i) = EP(ProblemFunction, DisplayFlag, false, PopSize, InitPopSize, RandSeed, MaxGen);
    MinCostAdapt(:, i) = EP(ProblemFunction, DisplayFlag, true, PopSize, InitPopSize, RandSeed, MaxGen);
end
MinCostNoAdapt = mean(MinCostNoAdapt, 2);
MinCostAdapt = mean(MinCostAdapt, 2);
figure
plot((0 : MaxGen), MinCostNoAdapt, 'r--', (0 : MaxGen), MinCostAdapt, 'b-')
xlabel('Generation')
ylabel('Minimum Cost')
legend('Without Variance Adaptation', 'With Variance Adaptation')