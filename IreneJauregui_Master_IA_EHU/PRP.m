function [Vnew, policy, states, Q1vector, Q2vector] = PRP(N, c1, c2)

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
tol = 10^(-5);

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
            else % then, Vberria(i)=Q2
                policy(i) = 2;
                if x2 == 0
                    Vnew(i) = 0;
                elseif y2 >= c2
                    Vnew(i) = 1;
                elseif N - x2 < c2 - y2
                    Vnew(i) = 0;
                elseif x2 < N && x2 > 0
                    Vnew(i) = 0;
                    for y = 0:x2
                        a = binopdf(y, x2, y2/x2) * binopdf(c2 - y2, N - x2, y/x2);
                        Vnew(i) = Vnew(i) + a;
                    end
                end
                continue;
            end
        end
        if x2 == N + 1
            if x1 == N
                Vnew(i) = (y1 >= c1);
                policy(i) = 0;
                continue;
            else % then, Vberria(i)=Q1
                policy(i) = 1;
                if x1==0
                    Vnew(i) = 0;
                elseif y1 >= c1
                    Vnew(i) = 1;
                elseif N - x1 < c1 - y1
                    Vnew(i) = 0;
                elseif x1 < N && x1 > 0
                    Vnew(i) = 0;
                    for y = 0:x1
                        a = binopdf(y, x1, y1/x1) * binopdf(c1 - y1, N - x1, y/x1);
                        Vnew(i) = Vnew(i) + a;
                    end
                end
                continue;
            end
        end

        % if N people have treatment 1 or 2
        if x1 == N && x2 == N
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
                if x2 == 0
                    Vnew(i) = 0;
                elseif y2 >= c2
                    Vnew(i) = 1;
                elseif N - x2 < c2 - y2
                    Vnew(i) = 0;
                elseif x2 < N && x2 > 0
                    Vnew(i) = 0;
                    for y = 0:x2
                        a = binopdf(y, x2, y2/x2) * binopdf(c2 - y2, N - x2, y/x2);
                        Vnew(i) = Vnew(i) + a;
                    end
                end
                continue;
            end
        end
        if x2 == N
            policy(i) = 1;
            if y2 >= c2
                Vnew(i) = 1;
                continue;
            else % then, Vnew(i)=Q1
                if x1==0
                    Vnew(i) = 0;
                elseif y1 >= c1
                    Vnew(i) = 1;
                elseif N - x1 < c1 - y1
                    Vnew(i) = 0;
                elseif x1 < N && x1 > 0
                    Vnew(i) = 0;
                    for y = 0:x1
                        a = binopdf(y, x1, y1/x1) * binopdf(c1 - y1, N - x1, y/x1);
                        Vnew(i) = Vnew(i) + a;
                    end
                end
                continue;
            end
        end
        
        % now x1 < N and x2 < N
        % 1ST TREATMENT
        if x1 == 0
            Q1 = 1;
        elseif y1 >= c1
            Q1 = 1;
        elseif N - x1 < c1 - y1
            Q1 = 0;
        elseif x1 < N && x1 > 0
            Q1 = 0;
            for y = 0:x1
                a = binopdf(y, x1, y1/x1) * binopdf(c1 - y1, N - x1, y/x1);
                Q1 = Q1 + a;
            end
        end
        

        % 2ND TREATMENT
        if x2 == 0
            Q2 = 1;
        elseif y2 >= c2
            Q2 = 1;
        elseif N - x2 < c2 - y2
            Q2 = 0;
        elseif x2 < N && x2 > 0
            Q2 = 0;
            for y = 0:x2
                a = binopdf(y, x2, y2/x2) * binopdf(c2 - y2, N - x2, y/x2);
                Q2 = Q2 + a;
            end
        end
                
        % SAVE THE BEST
        if abs(Q1-Q2) < tol
            Vnew(i) = Q1;
            policy(i) = 3;
        elseif Q1>Q2
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