function [T1,R1,T2,R2] = TR_aus_E(E)
% In dieser Funktion sollen die moeglichen euklidischen Transformationen
% aus der Essentiellen Matrix extrahiert werden
[U, Sigma, V] = svd(E);
RzPlusPi = [0, -1, 0;1, 0, 0;0, 0, 1];
RzMinusPi = [0, 1, 0;-1, 0, 0;0, 0, 1];
R1 = U * RzPlusPi' * V';
T1 = U * RzPlusPi * Sigma * U';
T1 = [-T1(2, 3), T1(1, 3), -T1(1, 2)]';

R2 = U * RzMinusPi' * V';
T2 = U * RzMinusPi * Sigma * U';
T2 = [-T2(2, 3), T2(1, 3), -T2(1, 2)]';
end