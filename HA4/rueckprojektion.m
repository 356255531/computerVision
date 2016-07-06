function repro_error = rueckprojektion(Korrespondenzen, P1, I2, T, R, K)
% Diese Funktion berechnet die projizierten Punkte in Kamera 2 und den
% mittleren Rueckprojektionsfehler
P1Homo = [P1;ones(1, size(P1, 2))];
P2Homo = [R, T;zeros(1, size(R, 2)), 1] * P1Homo;

x2Pro = K * [1 0 0 0;0 1 0 0;0 0 1 0] * P2Homo;

for i = 1 : size(x2Pro, 2)
	x2Pro(1:2, i) = x2Pro(1:2, i) / x2Pro(3, i);
end

x2BildKorrApp = x2Pro(1:2, :);
x2BildKorr = Korrespondenzen(3:4, :);
repro_error = mean(sqrt(sum((x2BildKorr - x2BildKorrApp) .^ 2)));

x2BildKorrApp = int16(x2BildKorrApp);
x2BildKorr = int16(x2BildKorr);
imshow(I2);
hold on;
scatter(x2BildKorr(1, :)', x2BildKorr(2, :)','+');
hold on;
scatter(x2BildKorrApp(1, :)', x2BildKorrApp(2, :)','*');
display(repro_error);
end