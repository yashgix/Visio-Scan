% paper=  corners = [173, 118; 509, 131;496, 607; 27, 530];
% paper2= corners= [295, 56; 526, 18; 582, 345; 347, 384];
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


maxR = round(sqrt(m^2 + n^2)); 
hLine = zeros(180, 2*maxR + 1);

for i = 1:m
    for j = 1:n
        if dilated_image(i, j) > 0
            for k = 1:180
                t = deg2rad(k - 1);
                R = round(j*cos(t) + i*sin(t)) + maxR + 1;
                hLine(k, R) = hLine(k, R) + 1;
            end
        end
    end
end


nSize = 4;
[row, col] = size(hLine);
hMax = false(row, col);

threshold = 0.53 * max(hLine(:)); 

for i = 2:row-1
    for j = 2:col-1
        neighborhood = hLine(i-1:i+1, j-1:j+1);
        max_val = max(neighborhood(:));
        if hLine(i, j) == max_val && max_val > threshold
            hMax(i, j) = true;
        end
    end
end

[maximaRows, maximaCols] = find(hMax);


intersection_points = [];

for i = 1:size(maximaRows, 1)-1
    for j = i+1:size(maximaRows, 1)
        theta1 = deg2rad(maximaRows(i) - 1);
        rho1 = maximaCols(i) - maxR - 1;
        theta2 = deg2rad(maximaRows(j) - 1);
        rho2 = maximaCols(j) - maxR - 1;
        x_intersect = (rho1*sin(theta2) - rho2*sin(theta1)) / sin(theta2 - theta1);
        y_intersect = (rho2*cos(theta1) - rho1*cos(theta2)) / sin(theta2 - theta1);
        intersection_points = [intersection_points; x_intersect, y_intersect];
    end
end


for k = 1:size(maximaRows, 1)
    R = maximaCols(k) - maxR - 1;
    t = deg2rad(maximaRows(k) - 1);
    x = 1:n;
    y = (R - x*cos(t)) / sin(t);
   
end

hold off;

detected_lines = [];

for i = 1:numel(maximaRows)
    R = maximaCols(i) - maxR - 1;
    t = maximaRows(i) - 1;

    detected_lines = [detected_lines; R, t];
end

%part 4

% Determine the four corner points of the paper

corner = [173, 118; 509, 131;496, 607; 27, 530 ];
ogCorners = corner;
rSize = [1100, 850]; 
rectImage = zeros(rSize(1), rSize(2), size(data, 3));
blank = [1, 1; rSize(2), 1; rSize(2), rSize(1); 1, rSize(1)];




I = [];
for i = 1:4
    X = ogCorners(i, 1);
    Y = ogCorners(i, 2);
    x = blank(i, 1);
    y = blank(i, 2);
    I = [I; -X, -Y, -1, 0, 0, 0, x*X, x*Y, x; 0, 0, 0, -X, -Y, -1, y*X, y*Y, y];
end

[U, S, V] = svd(I);
h = V(:, end);
homographyMat = reshape(h, [3, 3])';

for y = 1:rSize(1)
    for x = 1:rSize(2)
        ogLoc = homographyMat \ [x; y; 1];
        ogLoc = ogLoc(1:2) / ogLoc(3);
        if all(ogLoc >= 1) && all(ogLoc <= [size(data, 2); size(data, 1)])
            rectImage(y, x, :) = data(round(ogLoc(2)), round(ogLoc(1)), :);
        end
    end
end
figure;
imshow(data);
title('Original Image');

figure;
imshow(uint8(rectImage));
title('Rectified Image');