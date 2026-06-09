N = 15;
c1 = 3;
c2 = 3;
rho = 0.01;
epsilon = 0.1;
lambda = 0.99;
numepisodes = 10000000;
p1 = 0.1; %change also in Q-learning
p2 = 0.3;

% LOAD DATA
%data = load('N5_c1_p0.1_0.3.mat');
%d = data.datuak;
%states = d(:, 1:4);
%Qtable = d(:, 5:6);
%policy = d(:, 7);
%Q1vector1 = d(:, 8);
%Q2vector1 = d(:, 9);
%policy1 = d(:, 10);
%Q1vector2 = d(:, 11);
%Q2vector2 = d(:, 12);
%policy2 = d(:, 13);
%Q1vector3 = d(:, 14);
%Q2vector3 = d(:, 15);
%policy3 = d(:, 16);
%Q1vector4 = d(:, 17);
%Q2vector4 = d(:, 18);
%policy4 = d(:, 19);


% SAVE DATA
[states, Qtable, policy] = Qlearning(N, c1, c2, rho, epsilon, lambda, numepisodes);
[V1, policy1, states, Q1vector1, Q2vector1] = Value_iteration(N, lambda, c1, c2, p1, p2);
[V2, policy2, states, Q1vector2, Q2vector2] = Value_iteration_bayes(N, lambda, c1, c2, 1, 1);
[V3, policy3, states, Q1vector3, Q2vector3] = Oracle(N, c1, c2, p1, p2);
[V4, policy4, states, Q1vector4, Q2vector4] = PRP(N, c1, c2);
datuak = [states Qtable policy Q1vector1 Q2vector1 policy1 Q1vector2 Q2vector2 policy2 Q1vector3 Q2vector3 policy3 Q1vector4 Q2vector4 policy4];
save('N15_c3_p0.1_0.3.mat', 'datuak')