function [T,R, lambdas, P1] = rekonstruktion(T1,T2,R1,R2, Korrespondenzen, K)
% Funktion zur Bestimmung der korrekten euklidischen Transformation, der
% Tiefeninformation und der 3D Punkte der Merkmalspunkte in Bild 1
	[lambda1, P11] = get_lambda(Korrespondenzen, R1, T1, K);
	[lambda2, P12] = get_lambda(Korrespondenzen, R2, T2, K);
	if sum(find(lambda1 > 0)) > sum(find(lambda2 > 0))
		lambdas = lambda1;
		T = T1;
		R = R1;
        P1 = P11;
	else
		lambdas = lambda2;
		T = T2;
		R = R2;
        P1 = P12;
	end
end