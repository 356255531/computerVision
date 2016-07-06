%% get_lambda: get lambda from SVD
function [lambda, P] = get_lambda(Korrespondenzen, R, T, K)
	numSamples = size(Korrespondenzen, 2);

	x1 = [Korrespondenzen(1:2,:); ones(1,numSamples)];
	x2 = [Korrespondenzen(3:4,:); ones(1,numSamples)];
	x1 = K \ x1;
	x2 = K \ x2;
	x2Dach = get_dach(x2);

	diagEles = zeros(3, numSamples);
	for i = 1 : numSamples
		diagEles(:, i) = x2Dach(:, :, i) * R * x1(:, i);
	end

    Mleft = zeros(3 * numSamples, numSamples);
    for i = 1 : numSamples
    	Mleft(3 * i - 2:3 * i, i) = diagEles(:, i);
    end

    Mright = zeros(3 * numSamples, 1);
 	for i = 1 : numSamples
		Mright(i * 3 - 2:i * 3) = x2Dach(:, :, i) * T;
	end   

	M = [Mleft, Mright];

    [~, ~, V] = svd(M);
    d = V(:, end);
    lambda = d(1:end - 1);
	gamma = d(end);
	lambda = lambda / gamma;
	P = reshape(kron(lambda, [1,1,1]') .* x1(:), [3, size(x1, 2)]);
end
%% get_dach: get dach matrix expression
function [xDach] = get_dach(x)
	xDach = zeros(3, 3, size(x, 2));
	for i = 1 : size(x, 2)
		a = x(1, i);
		b = x(2, i);
		c = x(3, i);
		xDach(:, :, i) = [0, -c, b; c, 0, -a; -b, a, 0];
	end
end
