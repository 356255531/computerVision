function repro_error = rueckprojektion(Korrespondenzen, P1, I2, T, R, K)
% Diese Funktion berechnet die projizierten Punkte in Kamera 2 und den
% mittleren Rueckprojektionsfehler
P1Homo = [P1;ones(1, size(P1, 2))];
stdProMatrix = [1, 0, 0, 0;0, 1, 0, 0;0, 0, 1, 0];
xHomo = K * stdProMatrix * [R, T;zeros(1, size(R, 2)), 1] * P1Homo;
P2 = zeros(2, size(xHomo, 2));
for i = 1 : size(xHomo, 2)
	P2(:, i) = int16(xHomo(1:2, i) / xHomo(3, i));
end
hold on;
end