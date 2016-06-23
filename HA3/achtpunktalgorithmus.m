function [EF] = achtpunktalgorithmus(Korrespondenzen, K)
% Diese Funktion berechnet die Essentielle Matrix oder Fundamentalmatrix
% mittels 8-Punkt-Algorithmus, je nachdem, ob die Kalibrierungsmatrix 'K'
% vorliegt oder nicht
k = 8;

ifF = true;										% Ob K gegeben ist
if ~exist('K','var') || isempty(K)
    ifF = false;
end

if size(Korrespondenzen, 2) > k
	Korrespondenzen = Korrespondenzen(:, 1:k);	% When Korrespondenzen mehr als 8 ist, nur erste achte genommen sind.
end

x1 = [Korrespondenzen(1:2, :);ones(1, k)];
x1 = x1';
x2 = [Korrespondenzen(3:4, :);ones(1, k)];
x2 = x2';

A = zeros(k, 9);
for i  = 1:k
	A(i, :) = kron(x1(i, :), x2(i, :));			% Kron Produkt Berechnung
end

[~, ~, V] = svd(A);

EF = reshape(V(:,9), 3, 3)';

[U, D, V] = svd(EF);
D(3, 3) = 0;
EF = U * D * V';					% Essentielle Matrizen projektion

if ifF 											% When K gegeben ist, dann wird Fundamentalmatrizen berechnet.
	EF = K' * EF * K;
	[U,D,V] = svd(EF);
	sigma = (D(1, 1) + D(2, 2))/2;
	D(1,1) = sigma;
	D(2,2) = sigma;
	D(3,3) = 0;
	EF = U * D * V';
end
