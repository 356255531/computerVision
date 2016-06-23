function [Korrespondenzen_robust] = F_ransac(Korrespondenzen,varargin)
% Diese Funktion implementiert den RANSAC-Algorithmus zur Bestimmung von
% robusten Korrespondenzpunktpaaren

	%% Input parser
	P = inputParser;

	% Liste der optionalen Parameter
	% Die gesch??tzte Wahrscheilichkeit
	P.addOptional('epsilon', 0.65, @isnumeric)
	% Die gew??nschte Wahscheilichkeit
	P.addOptional('p', 0.99, @isnumeric);
	% Das Toleranzma??
	P.addOptional('tolerance',0.05, @isnumeric);

	% Lese den Input
	P.parse(varargin{:});

	epsilon 	= P.Results.epsilon;
	p 			= P.Results.p;
	tolerance 	= P.Results.tolerance;

	k = 8;
	iterZahl = log(1 - p) / log(1 - (1 - epsilon) ^ k);
	e3 = [0,-1,0;1,0,0;0,0,0];
	distSampson = @(x1, x2, F, e3) (x2' * F * x1) ^ 2 / (norm(e3 * F * x1) ^ 2 + norm(x2' * F * e3) ^ 2);

	minSum = inf;								% minimale Summe der Sampson-Distanzen
	for i = 1:iterZahl
		numKorr = size(Korrespondenzen, 2);		% Nummer des insgesamten Korrespondenzen
		idx = randperm(numKorr,k); 				% Random 8 Indexes der Korrespondenzen
		sample = Korrespondenzen(:,idx);		
		F = achtpunktalgorithmus(sample);		% Essentielle Matrizen Berechnung
		x1Pk = [sample(1:2, :);ones(1, k)];		% X1
		x2Pk = [sample(3:4, :);ones(1, k)];		% X2

		sumDist = 0;							% Summe der Sampson-Distanzen
		ifAccept = true;						% Ob alle 8 Sampson-Distanzen kleiner als Toleranz ist
		for j = 1:k
			x1 = x1Pk(:, j);
			x2 = x2Pk(:, j);
			dist = distSampson(x1, x2, F, e3);	% Berechnung derSampson-Distanzen

			if dist > tolerance
				ifAccept = false;
				break;
			else
				sumDist = sumDist + dist;
			end
		end

		if ifAccept && sumDist < minSum 		% Ob sich jetztige 8 Korrespondenzen die minimaleste Summe der Sampson-Distanzen verfuegt
			bestKorr = sample;
			minSum = sumDist;
		end

	end
	Korrespondenzen_robust	= bestKorr;

end