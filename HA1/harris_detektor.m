function  Merkmale = harris_detektor(Image, do_plot, segment_length, k, tau, N, min_dis) 
% In dieser Funktion soll der Harris-Detektor implementiert werden, der
% Merkmalspunkte aus dem Bild extrahiert

%Grauwertbild Beurteilung
if size(Image, 3) > 1
    fprintf('Es ist keine Grauwertbild!\n');
    fprintf('Das Bild wird ins Grauwertbild umgewandlt.');
    Image = rgb_to_gray(Image);
end

%Optionale Parameters Einstellung
if ~exist('do_plot','var') || isempty(do_plot)
    do_plot = 'False';
end
if ~exist('segment_length', 'var') || isempty(segment_length)
    segment_length = 50;
end
if ~exist('k', 'var') || isempty(k)
    k = 0.06;
end
if ~exist('tau','var') || isempty(tau)
    tau = 0.5;
end
if ~exist('N','var') || isempty(N)
    N = 3;
end
if ~exist('min_dis','var') || isempty(min_dis)
    min_dis = 10;
end

% Vorberechnung
sigma = 1;
sizeGaussianFilter = 3;
gaussianFilter = fspecial('gaussian',sizeGaussianFilter ^ 2, sigma); 
[gradientX, gradientY] = sobel_xy(Image);

gradientXY = conv2(gradientX .* gradientY, gaussianFilter, 'same');
gradientX2 = conv2(gradientX .^ 2, gaussianFilter, 'same');
gradientY2 = conv2(gradientY .^ 2, gaussianFilter,'same');

%H-Wert Berechnung G = [Gx2 Gxy;Gxy Gy2]
H = (gradientX2 .* gradientY2 - gradientXY .^ 2) - ...
    k * (gradientX2 + gradientY2) .^ 2;
H = ( 1000 / max(max(H))) * H;
H1 = H > tau;   % used as corner mark with value 0 or 1
H = H .* H1;

[m, n] = size(Image);

%Einteilung des Bildes in Segmente
for i = 1 : fix(m / segment_length)
    for j = 1 : fix(n / segment_length)
        mm = (i - 1) * segment_length;
        nn = (j - 1) * segment_length;

        H11 = H1(mm + 1 : mm + segment_length, nn + 1 : nn + segment_length);

        numEle = sum(sum(H11)); % Get the number of 1 element in H11

        if numEle < 2 % only one corner in segment then quit
            continue
        end


        HH = H(mm + 1 : mm + segment_length, nn + 1 : nn + segment_length); % get H value of segment

        [sortedHH, sortedIndex] = sort(HH(:), 'descend'); % do sort
        if numEle > N % reduce to N corner % Or when num of corner less then N, then take all
            sortedHH = sortedHH(1:N);
            sortedIndex = sortedIndex(1:N);
            H11 = HH > sortedHH(N);
        else
            sortedHH = sortedHH(1:numEle);
            sortedIndex = sortedIndex(1:numEle); 
        end

        c = fix((sortedIndex - 1) / segment_length) + 1; 
        r = mod(sortedIndex - 1, segment_length) + 1;

        segMerk = [r, c, sortedHH];
%%%%%%%%%%%% Perform minimal distance %%%%%%%%%%
        if size(segMerk, 1) > 1
            for x = 1 : size(segMerk,1) - 1
                if segMerk(x, 3) == 0
                    continue
                end

                for y = x + 1 : size(segMerk,1)
                    if segMerk(y, 3) == 0
                        continue
                    end

                    if norm(segMerk(x,1:2) - segMerk(y, 1:2)) < min_dis
                        segMerk(y, 3) = 0;
                        H11(segMerk(y, 1), segMerk(y, 2)) = 0;
                    end
                end
            end
        else % If N = 1, then take it
            H11 = zeros(segment_length, segment_length);
            H11(segMerk(1, 1), segMerk(1, 2)) = 1;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        H1(mm + 1 : mm + segment_length, nn + 1 : nn + segment_length) = H11;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Perform minimal distance twice%%%%%%%%%%%%%%%%%%
for i = 1 : fix(m / segment_length) - 1
    for j = 1 : fix(n / segment_length) - 1
        mm = fix((i - 0.5) * segment_length);
        nn = fix((j - 0.5) * segment_length);

        H11 = H1(mm + 1 : mm + segment_length, nn + 1 : nn + segment_length);

        numEle = sum(sum(H11));

        if numEle < 2
            continue
        end

        HH = H(mm + 1 : mm + segment_length, nn + 1 : nn + segment_length);

        [sortedHH, sortedIndex] = sort(HH(:), 'descend');

        sortedHH = sortedHH(1:numEle);
        sortedIndex = sortedIndex(1:numEle); 

        c = fix((sortedIndex - 1) / segment_length) + 1;
        r = mod(sortedIndex - 1, segment_length) + 1;

        segMerk = [r, c, sortedHH];

        for x = 1 : size(segMerk) - 1
            if segMerk(x, 3) == 0
                continue
            end

            for y = x + 1 : size(segMerk)
                if segMerk(y, 3) == 0
                    continue
                end

                if norm(segMerk(x,1:2) - segMerk(y, 1:2)) < min_dis
                    segMerk(y, 3) = 0;
                    H11(segMerk(y, 1), segMerk(y, 2)) = 0;
                end
            end
        end

        H1(mm + 1 : mm + segment_length, nn + 1 : nn + segment_length) = H11;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[r1, c1] = find(H1(11:end - 11, 11:end - 11)); % Take all corners
r1 = r1 + 5;
c1 = c1 + 5;
Merkmale = [r1, c1];


%Plot the Bild und markiert/ umrahmt den Pixel
if strcmp(do_plot, 'do_plot')
    sizeMerkmale = size(Merkmale, 1);
    halfDispLength = 5;
    for i = 1:sizeMerkmale
            Image(Merkmale(i, 1) : Merkmale(i, 1) + 2 * halfDispLength, ...
                Merkmale(i, 2)) = 255;
            Image(Merkmale(i, 1) : Merkmale(i, 1) + 2 * halfDispLength, ...
                Merkmale(i, 2) + 2 * halfDispLength) = 255;
            Image(Merkmale(i, 1), Merkmale(i, 2) : ...
                Merkmale(i, 2) + 2 * halfDispLength) = 255;
            Image(Merkmale(i, 1) + 2 * halfDispLength, Merkmale(i, 2) : ...
                Merkmale(i, 2) + 2 * halfDispLength) = 255;
    end
    imshow(uint8(Image));
end

end