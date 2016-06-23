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
	P.addOptional('tolerance', 0.2, @isnumeric);

	% Lese den Input
	P.parse(varargin{:});

	epsilon 	= P.Results.epsilon;
	p 			= P.Results.p;
	tolerance 	= P.Results.tolerance;

	k = 8;
	numKorr = size(Korrespondenzen, 2);			% Nummer des insgesamten Korrespondenzen
	iterZahl = log(1 - p) / log(1 - (1 - epsilon) ^ k);
	e3 = [0,-1,0;1,0,0;0,0,0];
	distSampson = @(x1, x2, F, e3) diag((x2' * F * x1) .^ 2) ./ (sum((e3 * F * x1) .^ 2)' + sum((x2' * F * e3)' .^ 2)');

	x1Pk = [Korrespondenzen(1:2, :);ones(1, numKorr)];		% X1
	x2Pk = [Korrespondenzen(3:4, :);ones(1, numKorr)];		% X2
	minSum = inf;								% Minimale Summe der Sampson-Distanzen
	maxConsNum = 0;								% Maximale nummer des Consensus Sets
	for i = 1:iterZahl
		idx = randperm(numKorr,k); 				% Random 8 Indexes der Korrespondenzen
		sample = Korrespondenzen(:,idx);		% Korrespondenzen benutzt in F Matrix Erzeugung	
		F = achtpunktalgorithmus(sample);		% Essentielle Matrizen Berechnung

		dist = distSampson(x1Pk, x2Pk, F, e3);	% Sampson Distanz vector

		idx = dist < tolerance;					% Index wo Sampson Distanz kleiner als tolerance ist
		cons = find(idx);						% Consensus Sets Index
		numCons = length(cons);
		sumDist = sum(dist(cons));				% Summe der Sampson-Distanzen

		if numCons > maxConsNum
			maxConsNum = numCons;
			minSum = sumDist;
			bestKorr = Korrespondenzen(:, cons);
        elseif numCons == maxConsNum
			if sumDist < minSum
				bestKorr = Korrespondenzen(:, cons);
				minSum = sumDist;
			end
		end
	end
	Korrespondenzen_robust	= bestKorr;

end