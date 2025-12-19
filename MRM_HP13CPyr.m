%% LV_ROI_thresholding_pipeline.m
%   1) Load a proton DICOM image (reference anatomy).
%   2) Interpolate proton image onto the 13C (pyruvate/metabolite) grid.
%   3) Draw LV Epicardial + Endocardial ROIs manually and derive myocardium.
%   4) Optionally extract a mid-myocardium ring via radial scaling.
%   5) Reconstruct / flip-angle-correct 13C images (pyruvate + metabolite).
%   6) Co-register 13C images to proton reference (translation model).
%   7) Compute sorted LV blood-pool signal distribution.
%   8) Sweep adaptive thresholds (alpha) and store thresholded LV voxel coords.
% Repository:
%   cor-utsw/hyperpolarized-13C
% Author:
%   Fatemeh Khashami
%   Dec 16, 2025
% Notes:
%   - This script is designed to be readable and modifiable for others.
%   - Several variables depend on the pipeline (see "USER INPUTS" section).
%   - Requires Image Processing Toolbox for interactive ROI drawing.
clc; clear; close all;
%%  USER INPUTS 
protonDicomPath = 'DICOM';  % replace with actual file path
% --- Geometric parameters (mm and matrix sizes) ---
pyruvate_FOV_mm        = 220;
pyruvate_recon_matrix  = 100;
proton_FOV_mm          = 350;
proton_recon_matrix    = 432;
% --- Mid-myocardium ring scaling (tune per subject) ---
scaling_factor_out = 0.10;   % outward expansion around mean radius
scaling_factor_in  = 0.10;   % inward contraction around mean radius
% --- LV endocardium ellipse radius (pixels on 100x100 interpolated proton) ---
newRadius_px = 9;            % adjust based on heart size
% --- Flip angles (degrees) ---
flipAngle_Pyr_deg = 3;       % typical: 3 or 5
flipAngle_Met_deg = 30;      % typical: 30 or 45
% --- Threshold sweep alpha for LV blood-pool detection ---
thresholds = 0.01:0.01:0.63;
%% Proton reference image used for ROI drawing and co-registration.
proton_Inj1 = dicomread(protonDicomPath);
protond     = double(proton_Inj1);
%%  INTERPOLATE PROTON -> 13C GRID 
% Proton is reconstructed at higher matrix size. We interpolate it onto the
% pyruvate reconstruction grid so ROIs drawn align with 13C image indexing.
pyruvate_grid_x = linspace(-pyruvate_FOV_mm/2, pyruvate_FOV_mm/2, pyruvate_recon_matrix);
pyruvate_grid_y = linspace(-pyruvate_FOV_mm/2, pyruvate_FOV_mm/2, pyruvate_recon_matrix);
proton_grid_x   = linspace(-proton_FOV_mm/2,   proton_FOV_mm/2,   proton_recon_matrix);
proton_grid_y   = linspace(-proton_FOV_mm/2,   proton_FOV_mm/2,   proton_recon_matrix);
[X,  Y ] = meshgrid(proton_grid_x,  proton_grid_y);        % size ~432x432
[Xq, Yq] = meshgrid(pyruvate_grid_x, pyruvate_grid_y);     % size ~100x100
interpolated_proton = interp2(X, Y, protond, Xq, Yq, 'linear', 0);
fixedImage = interpolated_proton;
%% DRAW ROIs (LV EPIC / ENDO)
% ROI #1: LV Epicardial contour (outer)
% ROI #2: LV Endocardial contour (inner)
% Myocardium should be: Epic minus Endo (i.e., outer ring).
figure('Name','Draw ROIs on interpolated proton (100x100)');
imagesc(fixedImage);
axis image; axis on;
try
    colormap turbo;
catch
    colormap parula;
end
title('Draw LV Epicardial ROI (outer) in RED, double-click to finalize');
binaryMasks = struct();
% --- ROI 1: Epicardium (outer) ---
h1 = imellipse(gca);
setColor(h1, 'r');
wait(h1);
binaryMask1 = createMask(h1);
binaryMasks.binaryMask1 = binaryMask1;
[row1, col1] = find(binaryMask1);
center_rc = [mean(row1), mean(col1)];  
% --- ROI 2: Endocardium (inner ---
title('Draw LV Endocardial ROI (inner) in BLUE, double-click to finalize');
h2 = imellipse(gca, [center_rc(2)-newRadius_px, center_rc(1)-newRadius_px, 2*newRadius_px, 2*newRadius_px]);
setColor(h2, 'b');
wait(h2);
binaryMask2 = createMask(h2);
binaryMasks.binaryMask2 = binaryMask2;
coordinates_mask1 = [row1, col1];          % epicardium coords
[row2, col2]       = find(binaryMask2);
coordinates_mask2  = [row2, col2];          % endocardium coords (LV blood pool ROI)
%% MYOCARDIUM 
myoMask = binaryMask1 & ~binaryMask2;
[rowM, colM]     = find(myoMask);
myo_coordinates  = [rowM, colM];
%%  MID-MYO RING (OPTIONAL)
% Create an annulus around the myocardium centroid using mean radius.
centroid_rc     = mean(myo_coordinates, 1);                 % [row, col]
distances       = sqrt(sum((myo_coordinates - centroid_rc).^2, 2));
mean_distance   = mean(distances);
inner_radius    = (1 - scaling_factor_in)  * mean_distance;
outer_radius    = (1 + scaling_factor_out) * mean_distance;
keep_inner      = distances > inner_radius;
keep_outer      = distances < outer_radius;
mid_myo_coordinates = myo_coordinates(keep_inner & keep_outer, :);
figure('Name','Myocardium + Mid-Myocardium Ring');
imagesc(fixedImage);
axis image; axis on;
try
    colormap turbo;
catch
    colormap parula;
end
hold on;
plot(myo_coordinates(:,2),      myo_coordinates(:,1),      'ob');  % myocardium (blue circles)
plot(mid_myo_coordinates(:,2),  mid_myo_coordinates(:,1),  '*r');  % mid-myo (red stars)
title('Myocardium (blue) and Mid-Myocardium Ring (red)');
hold off;
%% FLIP ANGLE CORRECTION FACTORS 
% Store sin(FA) factors separately to avoid overwriting degrees.
sinFA_Pyr = sin(deg2rad(flipAngle_Pyr_deg));
sinFA_Met = sin(deg2rad(flipAngle_Met_deg));
%% 13C DATA PLACEHOLDERS 
% IMPORTANT:
% The following variables must already exist in the workspace OR be loaded:
%   - img_sum_pyr            : 2D pyruvate summary image or frame (raw/combined)
%   - metabolic_PLB_inj1     : metabolite data array (e.g., [1, frame, x, y]%
% Example:
%   load('raw.mat'); load('workspace.mat'); etc.
Pyr_raw = img_sum_pyr(:,:);
extended_Pyr_raw = [Pyr_raw, zeros(size(Pyr_raw, 1), 1)];
flipped_Pyr_raw = flipud(extended_Pyr_raw);
metabolic_zerofillimg = double(zerofillimg(flipped_Pyr_raw, [100 100]));
metabolic_zerofillimg_Pyr_raw = abs(metabolic_zerofillimg).^2;
Pyr_FA =  metabolic_zerofillimg_Pyr_raw ./ sinFA_Pyr;
%% Bic for Frame 16, for example
B13Image(:,:) = metabolic_PLB_inj1(1, 16, :, :); %
extended_B13Image = [B13Image, zeros(size(B13Image, 1), 1)];
flip_carbonB_1 = flipud(extended_B13Image);
metabolic_zerofillimg_B_1 = double(zerofillimg(flip_carbonB_1, [100 100]));
metabolic_zerofillimgB2_1(:,:) = abs(metabolic_zerofillimg_B_1).^2;
Bic_FA =  metabolic_zerofillimgB2_1 ./ sinFA_Met;
%%  CO-REGISTRATION (TRANSLATION) 
% Estimate translation parameters to align 13C images to the proton reference.
% imgShifter is assumed to return [rowShift, colShift] or similar.
translationPlanP = imgShifter(metabolic_zerofillimg_Pyr_raw, fixedImage);
translationPlanB = imgShifter(metabolic_zerofillimgB2_1, fixedImage);
%%  SUM_Pyr PLACEHOLDER 
% SUM_Pyr is the sum over 44 frames of Pyr_FA{1,kkk1}.
matrixSize    = [pyruvate_recon_matrix, pyruvate_recon_matrix];
SUM_Pyr = zeros(matrixSize); 
for kkk1 = 1:44
    SUM_Pyr = SUM_Pyr + Pyr_FA{1, kkk1};
end
mask2_original_ROI = zeros(matrixSize);
mask2_original_ROI(sub2ind(matrixSize, coordinates_mask2(:,1), coordinates_mask2(:,2))) = 1;
%% sorted_signalPyr
signalPyr_Carbon_LV_ProtonROI = arrayfun(@(x, y) SUM_Pyr(y, x), coordinates_mask2(:, 2), coordinates_mask2(:, 1)); 
max_signalPyr_Carbon_LV_ProtonROI = max(signalPyr_Carbon_LV_ProtonROI);
min_signalPyr_Carbon_LV_ProtonROI = min(signalPyr_Carbon_LV_ProtonROI);
sorted_signalPyr = sort(signalPyr_Carbon_LV_ProtonROI, 'descend');
%% ADAPTIVE THRESHOLD SWEEP (α) 
% Threshold definition: alpha from 1% to 63%
%   threshold_E_1 = alpha * I_max 
%   cutoff        = threshold_E_1 + I_min
%   keep voxels where SUM_Pyr > cutoff
% Then restrict to LV mask2_original to keep only LV ROI voxels.
%% Cutoff formula for Blood pool ROI detection
for i = 1:length(thresholds)
   current_threshold = thresholds(i);
   threshold_E_1 = current_threshold * max_signalPyr_Carbon_LV_ProtonROI;
   max_min_LV = threshold_E_1 + min_signalPyr_Carbon_LV_ProtonROI;
   background_Thresh = SUM_Pyr > max_min_LV;
   threshold_E_1_all(i) = threshold_E_1;
   max_min_LV_all(i) = max_min_LV;
   background_Thresh_all{i} = background_Thresh;
   maskThresh_threshold = mask2_original_ROI & background_Thresh;
   [y_coords, x_coords] = find(maskThresh_threshold);
   coords_all{i} = [y_coords, x_coords];  % Exported ROI coordinates are applied to the pyruvate image for subsequent analysis.
end
%%
