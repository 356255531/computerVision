function [Korrespondenzen_robust] = F_ransac(Korrespondenzen,varargin)
% Diese Funktion implementiert den RANSAC-Algorithmus zur Bestimmung von
% robusten Korrespondenzpunktpaaren

	%% Input parser
	P = inputParser;

	% Liste der optionalen Parameter
	% Die gesch??tzte Wahrscheilichkeit
	P.addOptional('epsilon', 0.9, @isnumeric)
	% Die gew??nschte Wahscheilichkeit
	P.addOptional('p', 0.9, @isnumeric);
	% Das Toleranzma??
	P.addOptional('tolerance',0.9, @isnumeric);

	% Lese den Input
	P.parse(varargin{:});

	epsilon = P.Results.epsilon;
	p = P.Results.p;
	tolerance = P.Results.tolerance;

	k = 8;
	iterZahl = log(1 - p) / log(1 - (1 - epsilon) ^ k);
	e3 = [0,-1,0;1,0,0;0,0,0];
	distSampson = @(x1, x2, E, e) (x2' * E * x1) ^ 2 / (norm(e * E * x1) ^ 2 + norm(x2' * E * e) ^ 2);

	minSum = inf;
	for i = 1:iterZahl
		numKorr = size(Korrespondenzen, 2);
		idx = randperm(numKorr,k); 
		sample = Korrespondenzen(:,idx);
		E = achtpunktalgorithmus(sample);
		x1Pk = [sample(1:2, :);ones(1, k)];
		x2Pk = [sample(3:4, :);ones(1, k)];

		sumDist = 0;
		ifAccept = true;
		for j = 1:k
			x1 = x1Pk(:, j);
			x2 = x2Pk(:, j);
			dist = distSampson(x1, x2, E, e3);

			if dist > tolerance
				sumDist = sumDist + dist;
				ifAccept = false;
				break;
		end

		if ifAccept && sumDist < minSum
			bestKorr = sample;
			minSum = sumDist;
		end

	end

end