% POLITIKAK KONPARATZEKO, BALIO ITERAZIO ETA ORACLE METODOEN ARTEAN
% probabilitateak aldatzen direnean zein egoeretako politikak aldatzen dira?

N = 2;
lambda= 0.9999;
c1 = 2;
c2 = 1;

probabilitateak = 0:0.1:1;
bikoteak = [];
for i = 1:length(probabilitateak)
    for j = 1:length(probabilitateak)
        bikoteak = [bikoteak; probabilitateak(i), probabilitateak(j)];
    end
end
bikote_kopurua = size(bikoteak,1);

aldatutako_egoerak_Q1_Q2 = [];

for i = 1 : bikote_kopurua
    p1 = bikoteak(i,1);
    p2 = bikoteak(i,2);
    [~, politika, egoerak, Q1bektore, Q2bektore] = Value_iteration_treatments(N, lambda, c1, c2, p1, p2);
    [~, politika1, ~, Q1bektore1, Q2bektore1] = Oracle(N, c1, c2, p1, p2);
    egoera_kopurua = size(egoerak,1);
    for j = 1:egoera_kopurua
        if politika(j)~=politika1(j)
            aldatutako_egoerak_Q1_Q2 = [aldatutako_egoerak_Q1_Q2; egoerak(j,:) p1 p2 Q1bektore(j) Q2bektore(j) Q1bektore1(j) Q2bektore1(j)];
        end
    end
end