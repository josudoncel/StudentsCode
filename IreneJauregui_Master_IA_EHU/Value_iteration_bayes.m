function [Vnew, policy, states, Q1vector, Q2vector] = Value_iteration_bayes(N, lambda, c1, c2, alpha_prior, beta_prior)

% STATES IN ROWS
states = [];
for x1 = 0:N+1 % N+1, to take into account absorbing states
    for x2 = 0:N+1
        for y1 = 0:x1
            for y2 = 0:x2
                states = [states; x1 y1 x2 y2];
            end
        end
    end
end
number_states = size(states,1);

% INITIALIZATIONS
V = zeros(number_states,1);
Vnew = zeros(number_states,1);
policy = zeros(number_states,1);
Q1vector = zeros(number_states,1);
Q2vector = zeros(number_states,1);
dif = Inf;
tol = 10^(-8);

% MAIN LOOP
while dif > tol

    Vold = V;

    % analyze every state
    for i = 1:number_states
        x1 = states(i,1);
        y1 = states(i,2);
        x2 = states(i,3);
        y2 = states(i,4);

        p1 = (alpha_prior + y1)/(alpha_prior + beta_prior + x1);
        p2 = (alpha_prior + y2)/(alpha_prior + beta_prior + x2);
        
        % ABSORBING STATES
        % if we are already in a absorbing state
        absorbing_succes_index_1 = find(states(:,1)==x1+1 & states(:,2)==y1+1 & states(:,3)==x2 & states(:,4)==y2); %1 treatment success and 2 treatment in absorbing
        absorbing_failure_index_1 = find(states(:,1)==x1+1 & states(:,2)==y1 & states(:,3)==x2 & states(:,4)==y2); %1 treatment failure and 2 treatment in absorbing
        absorbing_succes_index_2 = find(states(:,1)==x1 & states(:,2)==y1 & states(:,3)==x2+1 & states(:,4)==y2+1); %2 treatment success and 1 treatment in absorbing
        absorbing_failure_index_2 = find(states(:,1)==x1 & states(:,2)==y1 & states(:,3)==x2+1 & states(:,4)==y2); %2 treatment failure and 1 treatment in absorbing
        if x1 == N + 1 && x2 == N + 1
            Vnew(i) = 0; 
            policy(i) = 0;
            continue;
        end
        if x1 == N + 1
            if x2 == N
                Q2vector(i) = (y2 >= c2) + lambda*(p2*Vold(absorbing_succes_index_2)+(1-p2)*Vold(absorbing_failure_index_2));
                Vnew(i) = Q2vector(i);
                policy(i) = 0;
                continue;
            else
                Q2vector(i) = lambda*(p2*Vold(absorbing_succes_index_2)+(1-p2)*Vold(absorbing_failure_index_2));
                Vnew(i) = Q2vector(i);
                policy(i) = 2;
                continue;
            end
        end
        if x2 == N + 1            
            if x1 == N
                Q1vector(i) = (y1 >= c1) + lambda*(p1*Vold(absorbing_succes_index_1)+(1-p1)*Vold(absorbing_failure_index_1));
                Vnew(i) = Q1vector(i);
                policy(i) = 0;
                continue;
            else
                Q1vector(i) = lambda*(p1*Vold(absorbing_succes_index_1)+(1-p1)*Vold(absorbing_failure_index_1));
                Vnew(i) = Q1vector(i);
                policy(i) = 1;
                continue;
            end
        end

        % if we jump now to a absorbing state
        to_absorbing_succes_index_1 = find(states(:,1)==x1+1 & states(:,2)==y1+1 & states(:,3)==x2+1 & states(:,4)==y2); %1 treatment success and 2 treatment to absorbing
        to_absorbing_failure_index_1 = find(states(:,1)==x1+1 & states(:,2)==y1 & states(:,3)==x2+1 & states(:,4)==y2); %1 treatment failure and 2 treatment to absorbing
        to_absorbing_succes_index_2 = find(states(:,1)==x1+1 & states(:,2)==y1 & states(:,3)==x2+1 & states(:,4)==y2+1); %2 treatment success and 1 treatment to absorbing
        to_absorbing_failure_index_2 = find(states(:,1)==x1+1 & states(:,2)==y1 & states(:,3)==x2+1 & states(:,4)==y2); %2 treatment failure and 1 treatment to absorbing
        if x1 == N && x2 == N % it wont happen, thus V=0
            Vnew(i) = 0; 
            policy(i) = 0;
            continue;
        end 
        if x1 == N
            Q2vector(i) = (y1 >= c1) + lambda*(p2*Vold(to_absorbing_succes_index_2)+(1-p2)*Vold(to_absorbing_failure_index_2));
            Vnew(i) = Q2vector(i);
            policy(i) = 2;
            continue;
        end
        if x2 == N
            Q1vector(i) = (y2 >= c2) + lambda*(p1*Vold(to_absorbing_succes_index_1)+(1-p1)*Vold(to_absorbing_failure_index_1));
            Vnew(i) = Q1vector(i);
            policy(i) = 1;
            continue;
        end

        % now x1 < N and x2 < N
        % TREATMENT 1
        succes_index_1 = find(states(:,1) == x1+1 & states(:,2) == y1+1 & states(:,3) == x2 & states(:,4) == y2);
        failure_index_1 = find(states(:,1) == x1+1 & states(:,2) == y1 & states(:,3) == x2 & states(:,4) == y2);

        Q1 = lambda*(p1*Vold(succes_index_1)+(1-p1)*Vold(failure_index_1));

        % TREATMENT 2
        succes_index_2 = find(states(:,1)==x1 & states(:,2)==y1 & states(:,3)==x2+1 & states(:,4)==y2+1);
        failure_index_2 = find(states(:,1)==x1 & states(:,2)==y1 & states(:,3)==x2+1 & states(:,4)==y2);

        Q2 = lambda*(p2*Vold(succes_index_2)+(1-p2)*Vold(failure_index_2));


        % SAVE THE BEST
        if abs(Q1 - Q2) < tol
            Vnew(i) = Q1;
            policy(i) = 3;
        elseif Q1 > Q2
            Vnew(i) = Q1;
            policy(i) = 1;
        elseif Q2 > Q1
            Vnew(i) = Q2;
            policy(i) = 2;
        end
        Q1vector(i) = Q1;
        Q2vector(i) = Q2;
    end

    dif = norm(Vnew-V, Inf);
    V = Vnew;
end
end