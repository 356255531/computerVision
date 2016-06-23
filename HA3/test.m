%  Gruppennummer:
%  Gruppenmitglieder:

%% Hausaufgabe 3
%  Bestimmung von robusten Korrespondenzpunktpaaren mittels RANSAC und
%  Berechnung der Essentiellen Matrix.

%  F?r die letztendliche Abgabe bitte die Kommentare in den folgenden Zeilen
%  enfernen und sicherstellen, dass alle optionalen Parameter ?ber den
%  entsprechenden Funktionsaufruf fun('var',value) modifiziert werden k?nnen.


%% Bilder laden
Image1 = imread('szeneL.jpg');
IGray1 = rgb_to_gray(Image1);

Image2 = imread('szeneR.jpg');
IGray2 = rgb_to_gray(Image2);

%% Harris-Merkmale berechnen
Merkmale1 = harris_detektor(IGray1,'segment_length',9,'k',0.05,'min_dist',80,'N',50,'do_plot',false);
Merkmale2 = harris_detektor(IGray2,'segment_length',9,'k',0.05,'min_dist',80,'N',50,'do_plot',false);

%% Korrespondenzsch?tzung
tic;
Korrespondenzen = punkt_korrespondenzen(IGray1,IGray2,Merkmale1,Merkmale2,'min_corr',0.92,'do_plot',false);
zeit_korrespondenzen = toc;
disp(['Es wurden ' num2str(size(Korrespondenzen,2)) ' Korrespondenzpunktpaare in ' num2str(zeit_korrespondenzen) 's gefunden.'])

Merkmale1 = Korrespondenzen(1:2, :)';
Merkmale2 = Korrespondenzen(3:4, :)';


E = estimateFundamentalMatrix(Merkmale1(1:16,:), Merkmale2(1:16, :));
[Merkmale2(1, :), 1] * E * [Merkmale1(1, :), 1]'
F = achtpunktalgorithmus(Korrespondenzen(:,1:8));
[Merkmale2(1, :), 1] * F * [Merkmale1(1, :), 1]'