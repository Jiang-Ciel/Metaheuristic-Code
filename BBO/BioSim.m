function BioSim(n)
% Run a biogeography simulation to calculate the probability of each species count.
% Modified November 26, 2013, to keep species counts in [0, n].
if ~exist('n', 'var') || isempty(n)
    n = 4;
end
mu = (0 : n) / n;
lambda = 1 - mu;
T = 10000 * (n+1);
S = round(n / 2);
SArray = zeros(1, T);
for i = 1 : T
    if rand < 0.5
        if rand < lambda(S+1)
            S = S + 1;
            S = min(S, n);
        end
    else
        if rand < mu(S+1)
            S = S - 1;
            S = max(S, 0);
        end
    end
    SArray(i) = S;
end
for i = 0 : n
    indices = find(SArray == i);
    disp(['Pr(S=', num2str(i), ') = ', num2str(length(indices) / length(SArray))]);
end