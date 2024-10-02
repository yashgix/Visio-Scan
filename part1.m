clear;
close all;
clc;

data = imread('paper.jpg');
x = rgb2gray(data);
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
Mag = smoothGXY;
smooth = atan2(smoothGy, smoothGx);

for i = 2:m-1
    for j = 2:n-1
        angle = rad2deg(smooth(i, j));
        if  (angle >= 7*pi/8 && angle < pi) ||  (angle >= 0 && angle < pi/8)
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
eImg(r2) = 0.03;

dilated_image = eImg;
structuring_element = ones(3);

for i = 2:m-1
    for j = 2:n-1
        if eImg(i, j) == 1
            dilated_image(i-1:i+1, j-1:j+1) = max(dilated_image(i-1:i+1, j-1:j+1), structuring_element);
        end
    end
end


figure;
imshow(data);
title('Original Image');


figure;
imshow(dilated_image);
title('Edge Image');














