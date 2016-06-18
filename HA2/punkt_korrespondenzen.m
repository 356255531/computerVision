function [Korrespondenzen] = punkt_korrespondenzen(I1,I2,Mpt1,Mpt2,inputs_args)
% In dieser Funktion sollen die extrahierten Merkmalspunkte aus einer
% Stereo-Aufnahme mittels NCC verglichen werden um Korrespondenzpunktpaare
% zu ermitteln.
[m, n] = size(I1);				% Size of image

Mpt1 = Mpt1';
Mpt2 = Mpt2';

numMPt1 = size(Mpt1, 1);					% Feature numbers in image1
numMPt2 = size(Mpt2, 1);					% Feature number in image2

window_length = inputs_args(1); 			% Window length
min_corr = inputs_args(2);					% Threashold
Korrespondenzen = [];

for i = 1 : numMPt1
	p1 = Mpt1(i, :);
    if p1(1) < fix(window_length / 2) + 1 || p1(1) > m - fix(window_length / 2) || ...
    	p1(2) < fix(window_length / 2) + 1 || p1(2) > n - fix(window_length / 2)
    	continue
    end

	a = I1(p1(1) - fix(window_length / 2) : ...
			p1(1) + fix(window_length / 2), ...
			p1(2) - fix(window_length / 2): ...
			p1(2) + fix(window_length / 2));

	a = double(a(:));
    mu = mean(a);
    a = a - mu;
	sigma = std(a);
	a = a / sigma;

	for j = 1:numMPt2
		p2 = Mpt2(j, :);
	    if p2(1) < fix(window_length / 2) + 1 || p2(1) > m - fix(window_length / 2) || ...
	    	p2(2) < fix(window_length / 2) + 1|| p2(2) > n - fix(window_length / 2)
	    	continue
	    end

		b = I2(p2(1) - fix(window_length / 2) : ...
				p2(1) + fix(window_length / 2), ...
				p2(2) - fix(window_length / 2): ...
				p2(2) + fix(window_length / 2));

		b = double(b(:));
		mu = mean(b);
        b = b - mu;
		sigma = std(b);
		b = b / sigma;

		korrK = a' * b / length(a) ;
		if korrK > min_corr

			Korrespondenzen = [Korrespondenzen, [p1' ; p2']];
		end
	end
end
showMatchedFeatures(I1,I2,Korrespondenzen(1:2,:)',Korrespondenzen(3:4,:)')
