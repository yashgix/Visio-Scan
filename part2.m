close all;
clear;
clc;

data = imread('paper.jpg'); 
x=im2gray(data);
x = double(x);

derX = [0, 0, 0;
        1/2, 0, -1/2;
        0, 0, 0];

derY = [0, 1/2, 0;
        0, 0, 0;
        0, -1/2, 0];

kernel = fspecial('gaussian', 3, 3);
M = conv2(x, kernel, 'same');

smoothGx = conv2(M, derX, 'same');
smoothGy = conv2(M, derY, 'same');
smoothGXY = sqrt(smoothGx .* smoothGx + smoothGy .* smoothGy);


[m, n] = size(smoothGXY);
sup = zeros(m, n);
Mag=smoothGXY;
smooth = atan2(smoothGy, smoothGx);

for i = 2:m-1
    for j = 2:n-1
        angle = rad2deg(smooth(i, j));
        if (angle >= 0 && angle < pi/8) || (angle >= 7*pi/8 && angle < pi)
            row1 = i;
            row2 = i;
            col1 = j + 1;
            col2 = j - 1;

        elseif (angle >= 3*pi/8 && angle < 5*pi/8)
            row1 = i - 1;
            row2 = i + 1;
            col1 = j;            
            col2 = j;

        elseif (angle >= pi/8 && angle < 3*pi/8)
            row1 = i + 1;
            row2 = i - 1;
            col1 = j - 1;
            col2 = j + 1;

        else
            row1 = i - 1;
            row2 = i + 1;
            col1 = j - 1;
            col2 = j + 1;
        end

        if (Mag(i, j) > Mag(row1, col1) && Mag(i, j) > Mag(row2, col2))
            sup(i, j) = Mag(i, j);
        end
    end
end

sup_max = max(sup(:));
th = 0.6 * sup_max;
tl = 0.2 * sup_max;
r1 = sup >= th;
r2 = (sup >= tl) & (sup < th);
eImg = zeros(size(sup));
eImg(r1) = 1;
eImg(r2) = 0.01;

[row, col] = find(eImg);
total = numel(row); 
maxR = sqrt(m^2 + n^2);  
hLine = zeros(round(2 * maxR), 180);
range = -90:180/179:90;

for i = 1:total
    x = col(i);
    y = row(i);
    for k = 1:180
        t = deg2rad(range(k)); 
        R = round((y * cos(t) - x * sin(t)) + maxR); 
        hLine(R, k) = hLine(R, k) + 1;
    end
end

hLine = hLine / max(hLine(:));  
figure;
imshow(hLine');
title('Hough Transform');



