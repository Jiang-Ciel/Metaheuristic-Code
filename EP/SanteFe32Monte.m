function SanteFe32Monte

nMonte = 1200;
DisplayFlag = false;
NumStates = 4 : 2 : 12;
FoodEaten = zeros(length(NumStates), nMonte);
Moves = zeros(length(NumStates), nMonte);
MaxFoodEaten = 0;
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    for StateNdx = 1 : length(NumStates)
        [FoodEaten(StateNdx, i), Moves(StateNdx, i), FSM] = SanteFe32(NumStates(StateNdx), DisplayFlag);
        if FoodEaten(StateNdx, i) > MaxFoodEaten
            MaxFoodEaten = FoodEaten(StateNdx, i);
            MinMoves = Moves(StateNdx, i);
            BestFSM = FSM;
        elseif (FoodEaten(StateNdx, i) == MaxFoodEaten) && (Moves(StateNdx, i) < MinMoves)
            MinMoves = Moves(StateNdx, i);
            BestFSM = FSM;
        end
    end
end
FoodEaten = mean(FoodEaten, 2);
disp('Average of food eaten for each FSM size:')
disp(FoodEaten)
disp(['Most food eaten = ', num2str(MaxFoodEaten)])
disp(['Fewest moves = ', num2str(MinMoves)])
disp('Best FSM:')
disp(BestFSM)