N = 15;
c1 = 3;
c2 = 3;
rho = 0.01;
epsilon = 0.1;
lambda = 0.99;
numepisodes = 10000000;
p1 = 0.1; %change also in Q-learning
p2 = 0.3;
S = 10000;


% LOAD DATA
data = load('N15_c3_p0.1_0.3.mat');
d = data.datuak;
states = d(:, 1:4);
policy = d(:, 7);
policy1 = d(:, 10);
policy2 = d(:, 13);
policy3 = d(:, 16);
policy4 = d(:, 19);




%Q-LEARNING
%[states, Qtable, policy] = Qlearning(N, c1, c2, rho, epsilon, lambda, numepisodes);

success = zeros(S, 2*N);
state_index = @(x) find(states(:,1)==x(1) & states(:,2)==x(2) & states(:,3)==x(3) & states(:,4)==x(4), 1);

for s = 1:S
    state = [0 0 0 0]; %initial state
    achieve1 = 0;
    achieve2 = 0;
    for i = 1:2*N

        if achieve1 && achieve2
            success(s,i) = 2;
            continue;
        end

        index = state_index(state);
        action = policy(index);

        if action == 3
            action = randi([1, 2]); %randomly
        end

        if action == 1
            a = (rand < p1);
            state = [state(1) + 1, state(2) + a, state(3), state(4)];
        end

        if action == 2
            a = (rand < p2);
            state = [state(1), state(2), state(3) + 1, state(4) + a];
        end

        if state(1) == N && state(2) >= c1
            achieve1 = 1;
        end

        if state(3) == N && state(4) >= c2
            achieve2 = 1;
        end
        success(s,i) = achieve1 + achieve2;
    end
end
success_q = mean(success, 1);




%OPTIMAL
%[~, policy1, states, ~, ~] = Value_iteration(N, lambda, c1, c2, p1, p2);

success = zeros(S, 2*N);
state_index = @(x) find(states(:,1)==x(1) & states(:,2)==x(2) & states(:,3)==x(3) & states(:,4)==x(4), 1);

for s = 1:S
    state = [0 0 0 0]; %initial state
    achieve1 = 0;
    achieve2 = 0;
    for i = 1:2*N

        if achieve1 && achieve2
            success(s,i) = 2;
            continue;
        end

        index = state_index(state);
        action = policy1(index);

        if action == 3
            action = randi([1, 2]); %randomly
        end

        if action == 1
            a = (rand < p1);
            state = [state(1) + 1, state(2) + a, state(3), state(4)];
        end

        if action == 2
            a = (rand < p2);
            state = [state(1), state(2), state(3) + 1, state(4) + a];
        end

        if state(1) == N && state(2) >= c1
            achieve1 = 1;
        end

        if state(3) == N && state(4) >= c2
            achieve2 = 1;
        end
        success(s,i) = achieve1 + achieve2;
    end
end
success_optimal = mean(success, 1);




%OPTIMAL BAYES
%[~, policy2, states, ~, ~] = Value_iteration_bayes(N, lambda, c1, c2, 1, 1);

success = zeros(S, 2*N);
state_index = @(x) find(states(:,1)==x(1) & states(:,2)==x(2) & states(:,3)==x(3) & states(:,4)==x(4), 1);

for s = 1:S
    state = [0 0 0 0]; %initial state
    achieve1 = 0;
    achieve2 = 0;
    for i = 1:2*N

        if achieve1 && achieve2
            success(s,i) = 2;
            continue;
        end

        index = state_index(state);
        action = policy2(index);

        if action == 3
            action = randi([1, 2]); %randomly
        end

        if action == 1
            a = (rand < p1);
            state = [state(1) + 1, state(2) + a, state(3), state(4)];
        end

        if action == 2
            a = (rand < p2);
            state = [state(1), state(2), state(3) + 1, state(4) + a];
        end

        if state(1) == N && state(2) >= c1
            achieve1 = 1;
        end

        if state(3) == N && state(4) >= c2
            achieve2 = 1;
        end
        success(s,i) = achieve1 + achieve2;
    end
end
success_bayes = mean(success, 1);



%ORACLE
%[Vberria3, policy3, states, Q1bektore3, Q2bektore3] = Oracle(N, c1, c2, p1, p2);


success = zeros(S, 2*N);
state_index = @(x) find(states(:,1)==x(1) & states(:,2)==x(2) & states(:,3)==x(3) & states(:,4)==x(4), 1);

for s = 1:S
    state = [0 0 0 0]; %initial state
    achieve1 = 0;
    achieve2 = 0;
    for i = 1:2*N

        if achieve1 && achieve2
            success(s,i) = 2;
            continue;
        end

        index = state_index(state);
        action = policy3(index);

        if action == 3
            action = randi([1, 2]); %randomly
        end

        if action == 1
            a = (rand < p1);
            state = [state(1) + 1, state(2) + a, state(3), state(4)];
        end

        if action == 2
            a = (rand < p2);
            state = [state(1), state(2), state(3) + 1, state(4) + a];
        end

        if state(1) == N && state(2) >= c1
            achieve1 = 1;
        end

        if state(3) == N && state(4) >= c2
            achieve2 = 1;
        end
        success(s,i) = achieve1 + achieve2;
    end
end
success_oracle = mean(success, 1);


%PRP
%[~, policy4, states, ~, ~] = PRP(N, c1, c2);

success = zeros(S, 2*N);
state_index = @(x) find(states(:,1)==x(1) & states(:,2)==x(2) & states(:,3)==x(3) & states(:,4)==x(4), 1);

for s = 1:S
    state = [0 0 0 0]; %initial state
    achieve1 = 0;
    achieve2 = 0;
    for i = 1:2*N

        if achieve1 && achieve2
            success(s,i) = 2;
            continue;
        end

        index = state_index(state);
        action = policy4(index);

        if action == 3
            action = randi([1, 2]); %randomly
        end

        if action == 1
            a = (rand < p1);
            state = [state(1) + 1, state(2) + a, state(3), state(4)];
        end

        if action == 2
            a = (rand < p2);
            state = [state(1), state(2), state(3) + 1, state(4) + a];
        end

        if state(1) == N && state(2) >= c1
            achieve1 = 1;
        end

        if state(3) == N && state(4) >= c2
            achieve2 = 1;
        end
        success(s,i) = achieve1 + achieve2;
    end
end
success_prp = mean(success, 1);


figure;
plot(0:2*N, [0, success_optimal], 'LineWidth', 2);
hold on
plot(0:2*N, [0, success_q], 'LineWidth', 2);
plot(0:2*N, [0, success_bayes], 'LineWidth', 2);
plot(0:2*N, [0, success_oracle], 'LineWidth', 2);
plot(0:2*N, [0, success_prp], 'LineWidth', 2);

legend('Optimal','Q-learning', 'Bayes', 'Oracle', 'PrP')
grid on;