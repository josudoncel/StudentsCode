function [Vnew, policy, states, Q1vector, Q2vector] = Oracle(N, c1, c2, p1, p2)

% STATES IN ROWS
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

% INITIALIZATIONS
V = ones(number_states,1);
Vnew = zeros(number_states,1);
policy = zeros(number_states,1);
Q1vector = zeros(number_states,1);
Q2vector = zeros(number_states,1);
dif = Inf;
tol = 10^(-8);

% MAIN LOOP
while dif > tol

    % analyze every state
    for i = 1:number_states 
        x1 = states(i,1);
        y1 = states(i,2);
        x2 = states(i,3);
        y2 = states(i,4);

        % ABSORBING STATES
        if x1 == N + 1 && x2 == N + 1
            Vnew(i) = 0; 
            policy(i) = 0;
            continue;
        end
        if x1 == N + 1
            if x2 == N
                Vnew(i) = (y2 >= c2);
                policy(i) = 0;
                continue;
            else
                Vnew(i) = 1 - binocdf(c2 - y2 - 1, N - x2, p2);
                policy(i) = 2;
                continue;
            end
        end
        if x2 == N + 1
            if x1 == N
                Vnew(i) = (y1 >= c1);
                policy(i) = 0;
                continue;
            else
                Vnew(i) = 1 - binocdf(c1 - y1 - 1, N - x1, p1);
                policy(i) = 1;
                continue;
            end
        end

        % if N people have treatment 1 or 2
        if x1 == N && x2 == N % it wont happen, thus V=0
            policy(i) = 0;
            Vnew(i) = 0;
            continue;
        end
        if x1 == N
            policy(i) = 2;
            if y1 >= c1
                Vnew(i) = 1;
                continue;
            else % then, Vnew(i)=Q2
                if y2 >= c2
                    Vnew(i) = 1;
                elseif N - x2 < c1 - y2
                    Vnew(i) = 0;
                else
                    Vnew(i) = 1 - binocdf(c2 - y2 - 1, N - x2, p2);
                end
                continue;
            end
        end
        if x2 == N
            policy(i) = 1;
            if y2 >= c2
                Vnew(i) = 1;
                continue;
            else % then, Vberria(i)=Q1
                if y1 >= c1
                    Vnew(i) = 1;
                elseif N - x1 < c1 - y1
                    Vnew(i) = 0;
                else
                    Vnew(i) = 1 - binocdf(c1 - y1 - 1, N - x1, p1);
                end
                continue;
            end
        end
        
        % now x1 < N and x2 < N
        % TREATMENT 1
        if y1 >= c1
            Q1 = 1;
        elseif N - x1 < c1 - y1
            Q1 = 0;
        else
            Q1 = 1 - binocdf(c1 - y1 - 1, N - x1, p1);
        end


        % TREATMENT 2      
        if y2 >= c2
            Q2 = 1;
        elseif N - x2 < c2 - y2
            Q2 = 0;
        else
            Q2 = 1 - binocdf(c2 - y2 - 1, N - x2, p2);
        end
        
        
        % SAVE THE BEST
        if abs(Q1-Q2) < tol
            Vnew(i) = Q1;
            policy(i) = 3;
        elseif Q1 > Q2
            Vnew(i) = Q1;
            policy(i) = 1;
        else
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