function [Vberria, politika, egoerak, Q1bektore, Q2bektore] = PRP(N, c1, c2)

% EGOERA POSIBLE GUZTIAK ERRENKADAKA IDATZI:
egoerak = [];
for x1 = 0:N+1
    for x2 = 0:N+1
        for y1 = 0:x1
            for y2 = 0:x2
                egoerak = [egoerak; x1 y1 x2 y2];
            end
        end
    end
end
egoera_kopurua = size(egoerak,1);

% ABIARAPENAK
V = ones(egoera_kopurua,1);
Vberria = zeros(egoera_kopurua,1);
politika = zeros(egoera_kopurua,1);
Q1bektore = zeros(egoera_kopurua,1);
Q2bektore = zeros(egoera_kopurua,1);
dif = Inf;
tol = 10^(-5);

% BUKLE NAGUSIA
while dif > tol

    % egoera bakoitza aztertu
    for i = 1:egoera_kopurua 
        x1 = egoerak(i,1);
        y1 = egoerak(i,2);
        x2 = egoerak(i,3);
        y2 = egoerak(i,4);
        
        % ABSORBING EGOERAK
        if x1 == N + 1 && x2 == N + 1
            Vberria(i) = 0; 
            politika(i) = 0;
            continue;
        end
        if x1 == N + 1
            if x2 == N
                Vberria(i) = (y2 >= c2);
                politika(i) = 0;
                continue;
            else % kasu honetan, Vberria(i)=Q2
                politika(i) = 2;
                if x2 == 0
                    Vberria(i) = 0;
                elseif y2 >= c2
                    Vberria(i) = 1;
                elseif N - x2 < c2 - y2
                    Vberria(i) = 0;
                elseif x2 < N && x2 > 0
                    Vberria(i) = 0;
                    for y = 0:x2
                        a = binopdf(y, x2, y2/x2) * binopdf(c2 - y2, N - x2, y/x2);
                        Vberria(i) = Vberria(i) + a;
                    end
                end
                continue;
            end
        end
        if x2 == N + 1
            if x1 == N
                Vberria(i) = (y1 >= c1);
                politika(i) = 0;
                continue;
            else % kasu honetan, Vberria(i)=Q1
                politika(i) = 1;
                if x1==0
                    Vberria(i) = 0;
                elseif y1 >= c1
                    Vberria(i) = 1;
                elseif N - x1 < c1 - y1
                    Vberria(i) = 0;
                elseif x1 < N && x1 > 0
                    Vberria(i) = 0;
                    for y = 0:x1
                        a = binopdf(y, x1, y1/x1) * binopdf(c1 - y1, N - x1, y/x1);
                        Vberria(i) = Vberria(i) + a;
                    end
                end
                continue;
            end
        end

        % N pertsonek 1 edo 2 tratamendua jaso badute 
        if x1 == N && x2 == N % berez hau ez da emango, orduan V=0
            politika(i) = 0;
            Vberria(i) = 0;
            continue;
        end
        if x1 == N
            politika(i) = 2;
            if y1 >= c1
                Vberria(i) = 1;
                continue;
            else % kasu honetan, Vberria(i)=Q2
                if x2 == 0
                    Vberria(i) = 0;
                elseif y2 >= c2
                    Vberria(i) = 1;
                elseif N - x2 < c2 - y2
                    Vberria(i) = 0;
                elseif x2 < N && x2 > 0
                    Vberria(i) = 0;
                    for y = 0:x2
                        a = binopdf(y, x2, y2/x2) * binopdf(c2 - y2, N - x2, y/x2);
                        Vberria(i) = Vberria(i) + a;
                    end
                end
                continue;
            end
        end
        if x2 == N
            politika(i) = 1;
            if y2 >= c2
                Vberria(i) = 1;
                continue;
            else % kasu honetan, Vberria(i)=Q1
                if x1==0
                    Vberria(i) = 0;
                elseif y1 >= c1
                    Vberria(i) = 1;
                elseif N - x1 < c1 - y1
                    Vberria(i) = 0;
                elseif x1 < N && x1 > 0
                    Vberria(i) = 0;
                    for y = 0:x1
                        a = binopdf(y, x1, y1/x1) * binopdf(c1 - y1, N - x1, y/x1);
                        Vberria(i) = Vberria(i) + a;
                    end
                end
                continue;
            end
        end
        
        % orain badakigu x1 < N eta x2 < N direla
        % 1. TRATAMENDUA
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
        

        % 2. TRATAMENDUA        
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
                
        % ONENA GORDE
        if abs(Q1-Q2) < tol
            Vberria(i) = Q1;
            politika(i) = 3;
        elseif Q1>Q2
            Vberria(i) = Q1;
            politika(i) = 1;
        else
            Vberria(i) = Q2;
            politika(i) = 2;
        end
        Q1bektore(i)=Q1;
        Q2bektore(i)=Q2;
    end
    dif = norm(Vberria-V, Inf);
    V=Vberria;
end
end