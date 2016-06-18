function [EF] = achtpunktalgorithmus(Korrespondenzen,K)
% Diese Funktion berechnet die Essentielle Matrix oder Fundamentalmatrix
% mittels 8-Punkt-Algorithmus, je nachdem, ob die Kalibrierungsmatrix 'K'
% vorliegt oder nicht

if ~exist('K','var') || isempty(K) || size(K) ~= [3, 3]
    ifF = false;
end

if size(Korrespondenzen, 2) > 8
	Korrespondenzen = Korrespondenzen(:, 8);
end

x1 = [Korrespondenzen(1:2, :);ones(1, 8)];
x1 = x1';
x2 = [Korrespondenzen(3:4, :);ones(1, 8)];
x2 = x2';

A = [];
for i  = 1:size(x1, 1)
	A = [A;kron(x1(i, :), x2(i, :))];
end

[U D V] = svd(A);

EF = reshape(V(:,9), 3, 3)';

[U D V] = svd(EF);

EF = U * diag([1, 1, 0]) * V';

if ifF
	inversedK = inv(K);
	EF = inversedK' * EF * inversedK;
end