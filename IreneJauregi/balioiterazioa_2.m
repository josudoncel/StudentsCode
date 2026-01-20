function [Vberria, politika, egoerak, Q1bektore, Q2bektore] = Value_iteration_treatments(N, lambda, c1, c2, p1, p2)
format long

% EGOERA POSIBLE GUZTIAK ERRENKADAKA IDATZI:
egoerak = [];
for x1 = 0:N+1 % N+1, absorbing egoerak ere kontuan hartzeko
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
V = zeros(egoera_kopurua,1);
Vberria = zeros(egoera_kopurua,1);
politika = zeros(egoera_kopurua,1);
Q1bektore = zeros(egoera_kopurua,1);
Q2bektore = zeros(egoera_kopurua,1);
dif = Inf;
tol = 10^(-8);

% BUKLE NAGUSIA
while dif > tol

    Vzaharra = V;

    % egoera bakoitza aztertu
    for i = 1:egoera_kopurua 
        x1 = egoerak(i,1);
        y1 = egoerak(i,2);
        x2 = egoerak(i,3);
        y2 = egoerak(i,4);
        
        % ABSORBING EGOERAK
        % aurretik absorbing egoeran bagaude
        absorbing_arrakasta_indizea_1 = find(egoerak(:,1)==x1+1 & egoerak(:,2)==y1+1 & egoerak(:,3)==x2 & egoerak(:,4)==y2); %1. tratamenduan arrakasta eta 2. tratamendua absorbing egoeran
        absorbing_porrot_indizea_1 = find(egoerak(:,1)==x1+1 & egoerak(:,2)==y1 & egoerak(:,3)==x2 & egoerak(:,4)==y2); %1. tratamenduan porrota eta 2. tratamendua absorbing egoeran
        absorbing_arrakasta_indizea_2 = find(egoerak(:,1)==x1 & egoerak(:,2)==y1 & egoerak(:,3)==x2+1 & egoerak(:,4)==y2+1); %2. tratamenduan arrakasta eta 1. tratamendua absorbing egoeran
        absorbing_porrot_indizea_2 = find(egoerak(:,1)==x1 & egoerak(:,2)==y1 & egoerak(:,3)==x2+1 & egoerak(:,4)==y2); %2. tratamenduan porrota eta 1. tratamendua absorbing egoeran
        if x1 == N + 1 && x2 == N + 1
            Vberria(i) = 0; 
            politika(i) = 0;
            continue;
        end
        if x1 == N + 1
            if x2 == N
                Vberria(i) = (y2 >= c2) + lambda*(p2*Vzaharra(absorbing_arrakasta_indizea_2)+(1-p2)*Vzaharra(absorbing_porrot_indizea_2));
                politika(i) = 0;
                continue;
            else
                Vberria(i) = lambda*(p2*Vzaharra(absorbing_arrakasta_indizea_2)+(1-p2)*Vzaharra(absorbing_porrot_indizea_2));
                politika(i) = 2;
                continue;
            end
        end
        if x2 == N + 1            
            if x1 == N
                Vberria(i) = (y1 >= c1) + lambda*(p1*Vzaharra(absorbing_arrakasta_indizea_1)+(1-p1)*Vzaharra(absorbing_porrot_indizea_1));
                politika(i) = 0;
                continue;
            else
                Vberria(i) = lambda*(p1*Vzaharra(absorbing_arrakasta_indizea_1)+(1-p1)*Vzaharra(absorbing_porrot_indizea_1));
                politika(i) = 1;
                continue;
            end
        end

        % orain sartu bagara absorbing egoerara
        absorbingera_arrakasta_indizea_1 = find(egoerak(:,1)==x1+1 & egoerak(:,2)==y1+1 & egoerak(:,3)==x2+1 & egoerak(:,4)==y2); %1. tratamenduan arrakasta eta 2. tratamendua absorbing egoerara
        absorbingera_porrot_indizea_1 = find(egoerak(:,1)==x1+1 & egoerak(:,2)==y1 & egoerak(:,3)==x2+1 & egoerak(:,4)==y2); %1. tratamenduan porrota eta 2. tratamendua absorbing egoerara
        absorbingera_arrakasta_indizea_2 = find(egoerak(:,1)==x1+1 & egoerak(:,2)==y1 & egoerak(:,3)==x2+1 & egoerak(:,4)==y2+1); %2. tratamenduan arrakasta eta 1. tratamendua absorbing egoerara
        absorbingera_porrot_indizea_2 = find(egoerak(:,1)==x1+1 & egoerak(:,2)==y1 & egoerak(:,3)==x2+1 & egoerak(:,4)==y2); %2. tratamenduan porrota eta 1. tratamendua absorbing egoerara
        if x1 == N && x2 == N % berez hau ez da emango, orduan V=0
            Vberria(i) = 0; 
            politika(i) = 0;
            continue;
        end 
        if x1 == N
            Vberria(i) = (y1 >= c1) + lambda*(p2*Vzaharra(absorbingera_arrakasta_indizea_2)+(1-p2)*Vzaharra(absorbingera_porrot_indizea_2));
            politika(i) = 2;
            continue;
        end
        if x2 == N
            Vberria(i) = (y2 >= c2) + lambda*(p1*Vzaharra(absorbingera_arrakasta_indizea_1)+(1-p1)*Vzaharra(absorbingera_porrot_indizea_1));
            politika(i) = 1;
            continue;
        end

        % orain badakigu x1 < N eta x2 < N direla
        % 1. TRATAMENDUA
        arrakasta_indizea_1 = find(egoerak(:,1) == x1+1 & egoerak(:,2) == y1+1 & egoerak(:,3) == x2 & egoerak(:,4) == y2);
        porrot_indizea_1 = find(egoerak(:,1) == x1+1 & egoerak(:,2) == y1 & egoerak(:,3) == x2 & egoerak(:,4) == y2);

        Q1 = lambda*(p1*Vzaharra(arrakasta_indizea_1)+(1-p1)*Vzaharra(porrot_indizea_1));

        % 2. TRATAMENDUA
        arrakasta_indizea_2 = find(egoerak(:,1)==x1 & egoerak(:,2)==y1 & egoerak(:,3)==x2+1 & egoerak(:,4)==y2+1);
        porrot_indizea_2 = find(egoerak(:,1)==x1 & egoerak(:,2)==y1 & egoerak(:,3)==x2+1 & egoerak(:,4)==y2);

        Q2=lambda*(p2*Vzaharra(arrakasta_indizea_2)+(1-p2)*Vzaharra(porrot_indizea_2));


        % ONENA GORDE
        if abs(Q1 - Q2) < tol
            Vberria(i) = Q1;
            politika(i) = 3;
        elseif Q1 > Q2
            Vberria(i) = Q1;
            politika(i) = 1;
        elseif Q2 > Q1
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