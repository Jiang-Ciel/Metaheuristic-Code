function PrimePredictionMonte(SizePenalty)
if ~exist('SizePenalty', 'var') || isempty(SizePenalty)
    SizePenalty = 0.5;
end
nMonte = 100;
GenLimit = 100;
MinCost = zeros(GenLimit+1, nMonte);
MinCostSizePenalty = zeros(GenLimit+1, nMonte);
NumStates = zeros(1, nMonte);
NumStatesSizePenalty = zeros(1, nMonte);
for i = 1 : nMonte
    disp(['Monte Carlo simulation ', num2str(i), '/', num2str(nMonte)])
    [NumStates(i), MinCost(:, i)] = PrimePrediction(GenLimit, 0, false);
    [NumStatesSizePenalty(i), MinCostSizePenalty(:, i)] = PrimePrediction(GenLimit, SizePenalty, false);
end
MinCost = mean(MinCost, 2);
MinCostSizePenalty = mean(MinCostSizePenalty, 2);
close all; SetPlotOptions
plot(0:GenLimit, MinCost, 'r--', 0:GenLimit, MinCostSizePenalty, 'b-')
xlabel('Generation')
ylabel('Minimum Cost')
legend('Without Size Penalty', 'With Size Penalty')
disp(['Average # of states: ', num2str(mean(NumStates)), ' without size penalty, ', ...
    num2str(mean(NumStatesSizePenalty)), ' with size penalty']);
disp(['Best performance = ', num2str(min(MinCost))]);