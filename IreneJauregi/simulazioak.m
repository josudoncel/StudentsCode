rng(0)
S = 10000; %simulazio kopurua

%OPTIMOA
[~, politika, egoerak, ~, ~] = Value_iteration_treatments(N, lambda, c1, c2, p1, p2);

arrakastak = zeros(S, 2*N);
egoera_indizea = @(x) find(egoerak(:,1)==x(1) & egoerak(:,2)==x(2) & egoerak(:,3)==x(3) & egoerak(:,4)==x(4), 1);

for s = 1:S
    egoera = [0 0 0 0]; %hasierako egoera
    lortu1 = 0;
    lortu2 = 0;
    for i = 1:2*N

        if lortu1 && lortu2
            arrakastak(s,i) = 2;
            continue;
        end

        indizea = egoera_indizea(egoera);
        akzioa = politika(indizea);

        if akzioa == 3
            akzioa = randi([1, 2]); %ausaz
        end

        if akzioa == 1
            a = (rand < p1);
            egoera = [egoera(1) + 1, egoera(2) + a, egoera(3), egoera(4)];
        end

        if akzioa == 2
            a = (rand < p2);
            egoera = [egoera(1), egoera(2), egoera(3) + 1, egoera(4) + a];
        end

        if egoera(1) == N && egoera(2) >= c1
            lortu1 = 1;
        end

        if egoera(3) == N && egoera(4) >= c2
            lortu2 = 1;
        end
        arrakastak(s,i) = lortu1 + lortu2;
    end
end
arrakastak_bbt_optimoa = mean(arrakastak, 1);



%ORACLE
[~, politika1, ~, ~, ~] = Oracle(N, c1, c2, p1, p2);

arrakastak = zeros(S, 2*N);
egoera_indizea = @(x) find(egoerak(:,1)==x(1) & egoerak(:,2)==x(2) & egoerak(:,3)==x(3) & egoerak(:,4)==x(4), 1);

for s = 1:S
    egoera = [0 0 0 0]; %hasierako egoera
    lortu1 = 0;
    lortu2 = 0;
    for i = 1:2*N

        if lortu1 && lortu2
            arrakastak(s,i) = 2;
            continue;
        end

        indizea = egoera_indizea(egoera);
        akzioa = politika1(indizea);

        if akzioa == 3
            akzioa = randi([1, 2]); %ausaz
        end

        if akzioa == 1
            a = (rand < p1);
            egoera = [egoera(1) + 1, egoera(2) + a, egoera(3), egoera(4)];
        end

        if akzioa == 2
            a = (rand < p2);
            egoera = [egoera(1), egoera(2), egoera(3) + 1, egoera(4) + a];
        end

        if egoera(1) == N && egoera(2) >= c1
            lortu1 = 1;
        end

        if egoera(3) == N && egoera(4) >= c2
            lortu2 = 1;
        end
        arrakastak(s,i) = lortu1 + lortu2;
    end
end
arrakastak_bbt_oracle = mean(arrakastak, 1);




%PRP
[~, politika2, ~, ~, ~] = PRP(N, c1, c2);

arrakastak = zeros(S, 2*N);
egoera_indizea = @(x) find(egoerak(:,1)==x(1) & egoerak(:,2)==x(2) & egoerak(:,3)==x(3) & egoerak(:,4)==x(4), 1);

for s = 1:S
    egoera = [0 0 0 0]; %hasierako egoera
    lortu1 = 0;
    lortu2 = 0;
    for i = 1:2*N

        if lortu1 && lortu2
            arrakastak(s,i) = 2;
            continue;
        end

        indizea = egoera_indizea(egoera);
        akzioa = politika2(indizea);

        if akzioa == 3
            akzioa = randi([1, 2]); %ausaz
        end

        if akzioa == 1
            a = (rand < p1);
            egoera = [egoera(1) + 1, egoera(2) + a, egoera(3), egoera(4)];
        end

        if akzioa == 2
            a = (rand < p2);
            egoera = [egoera(1), egoera(2), egoera(3) + 1, egoera(4) + a];
        end

        if egoera(1) == N && egoera(2) >= c1
            lortu1 = 1;
        end

        if egoera(3) == N && egoera(4) >= c2
            lortu2 = 1;
        end
        arrakastak(s,i) = lortu1 + lortu2;
    end
end
arrakastak_bbt_prp = mean(arrakastak, 1);

figure;
plot(0:2*N, [0, arrakastak_bbt_optimoa], 'LineWidth', 2);
hold on
plot(0:2*N, [0, arrakastak_bbt_oracle], 'LineWidth', 2);
plot(0:2*N, [0, arrakastak_bbt_prp], 'LineWidth', 2);
legend('Optimoa','Oracle','PrP')
grid on;