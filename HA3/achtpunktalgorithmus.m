function [EF] = achtpunktalgorithmus(Korrespondenzen, K)
% Diese Funktion berechnet die Essentielle Matrix oder Fundamentalmatrix
% mittels 8-Punkt-Algorithmus, je nachdem, ob die Kalibrierungsmatrix 'K'
% vorliegt oder nicht
numSamples = size(Korrespondenzen, 2);			% Nummer der Korrespondenzen

x1 = [Korrespondenzen(1:2, :);ones(1, numSamples)];
x1 = x1';
x2 = [Korrespondenzen(3:4, :);ones(1, numSamples)];
x2 = x2';

A = zeros(numSamples, 9);
for i  = 1:numSamples
	A(i, :) = kron(x1(i, :), x2(i, :));			% Kron Produkt Berechnung
end

[~, ~, V] = svd(A);
EF = reshape(V(:, 9), 3, 3);		

[U, D, V] = svd(EF);
D(3, 3) = 0;
EF = U * D * V';					% Essentielle Matrizen projektion

if exist('K','var') && ~isempty(K)		% When K gegeben ist, dann wird Fundamentalmatrizen berechnet.
	EF = K' * EF * K;
	[U,D,V] = svd(EF);
	sigma = (D(1, 1) + D(2, 2))/2;
	EF = U * diag([sigma, sigma, 0]) * V';
end