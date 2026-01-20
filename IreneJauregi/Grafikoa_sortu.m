function grafikoa = Grafikoa_sortu(politika_matrizea, bikoteak, egoerak, aldatutako_egoera)

% GRAFIKOKO POLTIKAK
% 0 --> .
% 1 --> x
% 2 --> +
% 3 --> *

egoera_indizea = find(ismember(egoerak, aldatutako_egoera, 'rows'), 1); % ismemberrek egoerak matrizeko ilera bakoitza aldatutako_egoera den ikusiko du 

% EGOERETAKO POLITIKAK
politikak = politika_matrizea(:, egoera_indizea);
p1 = bikoteak(:,1);
p2 = bikoteak(:,2);

% AKZIOEN INDIZEAK
indizea0 = find(politikak == 0);
indizea1 = find(politikak == 1);
indizea2 = find(politikak == 2);
indizea3 = find(politikak == 3);

% GRAFIKOA SORTU 
grafikoa = figure;
hold on;
g1 = plot(p1(indizea0), p2(indizea0), 'k.', 'MarkerSize', 12);
g2 = plot(p1(indizea1), p2(indizea1), 'kx', 'MarkerSize', 8, 'LineWidth', 1.2);
g3 = plot(p1(indizea2), p2(indizea2), 'k+', 'MarkerSize', 8, 'LineWidth', 1.2);
g4 = plot(p1(indizea3), p2(indizea3), 'k*', 'MarkerSize', 8, 'LineWidth', 1.2);
xlabel('p1');
ylabel('p2');
title(sprintf('[%d %d %d %d] egoerarako politiken grafikoa ', aldatutako_egoera(1), aldatutako_egoera(2), aldatutako_egoera(3), aldatutako_egoera(4)));
xlim([-0.05 1.05]);
ylim([-0.05 1.05]);
axis square;
grid on;
legend([g1, g2, g3, g4], {'0 --> .', '1 --> x', '2 --> +', '3 --> *'});
hold off;

xlabel('p1');
ylabel('p2');
title(sprintf('[%d %d %d %d] egoerarako politiken grafikoa ', aldatutako_egoera(1), aldatutako_egoera(2), aldatutako_egoera(3), aldatutako_egoera(4)));
xlim([-0.05 1.05]);
ylim([-0.05 1.05]);
axis square;
grid on;

% LEGENDA
legend_posibleak = {'0', '1', '2', '3'};
a = [];
legenda = {};
if ~isempty(g1)
    a = [a, g1];
    legenda = [legenda, legend_posibleak{1}];
end
if ~isempty(g2)
    a = [a, g2];
    legenda = [legenda, legend_posibleak{2}];
end
if ~isempty(g3)
    a = [a, g3];
    legenda = [legenda, legend_posibleak{3}];
end
if ~isempty(g4)
    a = [a, g4];
    legenda = [legenda, legend_posibleak{4}];
end
legend(a, legenda, 'Location', 'northeastoutside');

hold off;
end