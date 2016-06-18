function [Gray_image] = rgb_to_gray(Image)
% Diese Funktion soll ein RGB-Bild in ein Graustufenbild umwandeln. Falls
% das Bild bereits in Graustufen vorliegt, soll es direkt zur?ckgegeben werden.
if size(Image, 3) ~= 1 
    Image = double(Image);
	weight = [0.299 0.587 0.114];
	weightedR = Image(:,:,1) * weight(1);
	weightedG = Image(:,:,2) * weight(2);
	weightedB = Image(:,:,3) * weight(3);
	Gray_image = weightedR + weightedG + weightedB;
else
    Gray_image = double(Image);
end