function [states, Qtable, policy] = Qlearning(N, c1, c2, rho, epsilon, lambda, numepisodes)

p1 = 0.1;
p2 = 0.3;

states = [];
for x1 = 0:N+1
    for x2 = 0:N+1
        for y1 = 0:x1
            for y2 = 0:x2
                states = [states; x1 y1 x2 y2];
            end
        end
    end
end
number_states = size(states,1);
no_absorbing_state = find(states(:,1) < N & states(:,3) < N);

Qtable = zeros(number_states, 2);
print = false;

for ep = 1:numepisodes
    i = no_absorbing_state(randi(length(no_absorbing_state)));
    s = states(i, :);
    a0 = [];
    state_index0 = [];
    firststep = true; %in the first step, it doesnt update, because there is no previous state
    s0 = [];

    if print
        ep
        s
    end

    while s(1) < N+1 || s(3) < N+1
        x1 = s(1);
        y1 = s(2);
        x2 = s(3);
        y2 = s(4);
        state_index = find(states(:,1)==x1 & states(:,2)==y1 & states(:,3)==x2 & states(:,4)==y2, 1);

        if rand < epsilon
            if x1 < N && x2 < N
                a = randi([1 2]);
            elseif x1 < N
                a = 1;
            elseif x2 < N
                a = 2;
            end
        else
            if x1 < N && x2 < N
                Q1 = Qtable(state_index, 1);
                Q2 = Qtable(state_index, 2);
                if Q1 > Q2
                    a = 1;
                elseif Q2 > Q1
                    a = 2;
                else
                    a = randi([1 2]);
                end
            elseif x1 < N
                a = 1;
            elseif x2 < N
                a = 2;
            end
        end

        x1new = x1;
        y1new = y1;
        x2new = x2;
        y2new = y2;

        if a == 1
            x1new = x1new + 1;
            if rand < p1
                y1new = y1new + 1;
            end
        elseif a == 2
            x2new = x2new + 1;
            if rand < p2
                y2new = y2new + 1;
            end
        end

        if a == 2 && x1 == N
            x1new = N + 1;
        end
        if a == 1 && x2 == N
            x2new = N + 1;
        end

        snew = [x1new, y1new, x2new, y2new];

        if print
            a
            snew
        end

        if firststep
            firststep = false;
        else
            x1_0 = s0(1);
            y1_0 = s0(2);
            x2_0 = s0(3);
            y2_0 = s0(4);
            r = (y1_0 >= c1 && x1_0 == N) || (y2_0 >= c2 && x2_0 == N);
            Qtable(state_index0, a0) = Qtable(state_index0, a0) + rho * (r + lambda * max(Qtable(state_index, :)) - Qtable(state_index0, a0));
            if print
                s0
                r
                Qtable
            end
        end

        state_index0 = state_index;
        a0 = a;
        s0 = s;
        s = snew;

        if s(1) == N+1 && s(3) == N+1
            break;
        end
    end

    % Last update
    x1_0 = s0(1);
    y1_0 = s0(2);
    x2_0 = s0(3);
    y2_0 = s0(4);
    r = (y1_0 >= c1 && x1_0 == N) || (y2_0 >= c2 && x2_0 == N);
    x1 = s(1); 
    y1 = s(2); 
    x2 = s(3); 
    y2 = s(4);
    state_index_fin = find(states(:,1)==x1 & states(:,2)==y1 & states(:,3)==x2 & states(:,4)==y2, 1);
    Qtable(state_index0, a0) = Qtable(state_index0, a0) + rho * (r + lambda * max(Qtable(state_index_fin, :)) - Qtable(state_index0, a0));

    if print
        r
        Qtable
    end

end

% Equalize equivalent states
for i = 1:number_states
    x1 = states(i,1);
    y1 = states(i,2);
    x2 = states(i,3);
    y2 = states(i,4);

    % (0,0,N+1,y2) -> (0,0,N,min(y2,N)) - (y2>=c2)
    if x1 == 0 && x2 == N+1
        y2_eq = min(y2, N);
        i_eq = find(states(:,1)==0 & states(:,2)==0 & states(:,3)==N & states(:,4)==y2_eq, 1);
        if ~isempty(i_eq)
            Qtable(i,1) = Qtable(i_eq,1) - (y2_eq >= c2);
        end
    end

    % (N+1,y1,0,0) -> (N,min(y1,N),0,0) - (y1>=c1)
    if x2 == 0 && x1 == N+1
        y1_eq = min(y1, N);
        i_eq = find(states(:,1)==N & states(:,2)==y1_eq & states(:,3)==0 & states(:,4)==0, 1);
        if ~isempty(i_eq)
            Qtable(i,2) = Qtable(i_eq,2) - (y1_eq >= c1);
        end
    end

    % If x2 = N+1, (x1, y1, N+1, y2) -> (x1, y1, N+1, 0)
    if x2 == N+1
        i2 = find(states(:,1)==x1 & states(:,2)==y1 & states(:,3)==N+1 & states(:,4)==0, 1);
        Qtable(i, :) = Qtable(i2, :);
    end
    
    % If x1 = N+1, (N+1, y1, x2, y2) -> (N+1, 0, x2, y2)
    if x1 == N+1
        i1 = find(states(:,1)==N+1 & states(:,2)==0 & states(:,3)==x2 & states(:,4)==y2, 1);
        Qtable(i, :) = Qtable(i1, :);
    end
end

% SAVE POLICY
policy = zeros(number_states,1);
for i = 1:number_states
    x1 = states(i,1);
    x2 = states(i,3);
    if x1 < N || x2 < N
        if abs(Qtable(i,2) - Qtable(i,1)) < 10^(-8) && Qtable(i,2) ~= 0
        elseif Qtable(i,2) < Qtable(i,1)
            policy(i) = 1;
        elseif Qtable(i,2) > Qtable(i,1)
            policy(i) = 2;
        end
    end
end
end