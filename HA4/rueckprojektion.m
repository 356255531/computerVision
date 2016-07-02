function repro_error = rueckprojektion(Korrespondenzen, P1, I2, T, R, K)
% Diese Funktion berechnet die projizierten Punkte in Kamera 2 und den
% mittleren Rueckprojektionsfehler
P1Homo = [P1;ones(1, size(P1, 2))];
P2Homo = [R, T;zeros(1, size(R, 2)), 1] * P1Homo;

newKorrespondenzen = [Korrespondenzen(3:4, :);Korrespondenzen(1:2, :)];
E = achtpunktalgorithmus(newKorrespondenzen, K);
[T1,R1,T2,R2] = TR_aus_E(E);
[~, ~, lambda, ~] = rekonstruktion(T1,T2,R1,R2,newKorrespondenzen,K);
P2App = P2Homo(1:3, :);
P2Pro = reshape(P2App(:) ./ kron(lambda, [1,1,1]'), [3, size(P2App, 2)]);
x2Pro = K * P2Pro;
hold on;
end