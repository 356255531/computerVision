function [Fx,Fy] = sobel_xy(Image)
% In dieser Funktion soll das Sobel-Filter implementiert werden, welches
% ein Graustufenbild einliest und den Bildgradienten in x- sowie in
% y-Richtung zur?ckgibt.
if size(Image, 3) > 1
    Image = rgb_to_gray(Image);
end

Image = double(Image);
x = [-1 0 1;-2 0 2;-1 0 1];
y = x';

Fx = conv2(Image, x, 'same');
Fy = conv2(Image, y, 'same');