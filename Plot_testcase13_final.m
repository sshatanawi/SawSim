
% MATLAB Script: plot_gw_results.m
clear variables;

% 1) Read the CSV Data and Build the H Array
data = readmatrix('GW_results.csv');  % columns: [step, i, j, k, head]

desired_step  = 2;    % e.g., time step = 1
desired_layer = 2;    % e.g., layer = 1
dx = 500;             % cell width (ft)
dy = 500;             

% Filter for the desired time step and layer
rows = (data(:,1) == desired_step) & (data(:,4) == desired_layer);
selectedData = data(rows, :);

% Determine grid dimensions
NCOL = max(selectedData(:,2));
NROW = max(selectedData(:,3));

% Build H array (rows: 1..NROW, columns: 1..NCOL)
H = nan(NROW, NCOL);
noFlowMarker = 999.999;
tol = 1e-4;

for r = 1:size(selectedData,1)
    i_col   = selectedData(r,2);  % column index
    j_row   = selectedData(r,3);  % row index
    headVal = selectedData(r,5);
    if abs(headVal - noFlowMarker) < tol
        headVal = NaN;
    end
    H(j_row, i_col) = headVal;
end

% 2) Define Common Parameters for Both Plots

% For pcolor, define cell edges (NCOL+1 and NROW+1 edges)
x_edges = 0:dx:(NCOL*dx);
y_edges = 0:dy:(NROW*dy);

% For contour plots, define cell centers
x_centers = ((1:NCOL) - 0.5) * dx;
y_centers = ((1:NROW) - 0.5) * dy;
[Xc, Yc] = meshgrid(x_centers, y_centers);

% Define well locations (cell centers)
x_well1 = (4 - 0.5) * dx;  % well at cell (4,3)
y_well1 = (3 - 0.5) * dy;
x_well2 = (4 - 0.5) * dx;  % well at cell (4,8)
y_well2 = (8 - 0.5) * dy;

% Define river parameters:
% River along column 14 for rows 2 to 9
river_col = 14;
river_rows = 2:9;
x_river_center = (river_col - 0.5) * dx;
y_river = (river_rows - 0.5) * dy;
% Define patch vertices for the river rectangle.
% The river patch spans from (x_river_center - dx/2) to (x_river_center + dx/2)
% and from the lower edge of row 2 to the upper edge of row 9.
river_x = [x_river_center - dx/2, x_river_center + dx/2, ...
           x_river_center + dx/2, x_river_center - dx/2];
river_y = [min(y_river) - dy/2, min(y_river) - dy/2, ...
           max(y_river) + dy/2, max(y_river) + dy/2];

% 3) Plot 1: pcolor with Contour Overlay

% Pad H so that pcolor shows all data values
H_pad = nan(NROW+1, NCOL+1);
H_pad(1:end-1, 1:end-1) = H;
H_pad(end, 1:end-1) = H(end,:);
H_pad(1:end-1, end) = H(:, end);
H_pad(end, end) = H(end,end);

figure;
pcolor(x_edges, y_edges, H_pad);
shading flat;           % Remove grid lines between patches
colormap('summer');
cb = colorbar;
cb.Label.String = 'Groundwater Head (ft)';
cb.FontSize = 12;
cb.FontName = 'Times New Roman';
cb.Label.FontSize = 12;
cb.Label.FontName = 'Times New Roman';
xlabel('x-coordinate (ft)');
ylabel('y-coordinate (ft)');
%title('A. Upper aquifer under natural conditions');
% title('B. Lower aquifer under natural conditions');
% title('C. Upper aquifer under pumping conditions ');
 title('D. Lower Aquifer under pumping conditions ');
axis equal;
set(gca, 'YDir', 'normal');  % y increases upward
grid on;
hold on;

% Overlay contour lines using the original H data
contour(Xc, Yc, H, 'LineColor', 'k', 'LineWidth', 1.5, 'ShowText', 'on');

% Overlay wells using scatter with NO transparency
% hWell1 = scatter(x_well1, y_well1, 200, 'r', 'filled', 'DisplayName', 'Well');
% hWell1.MarkerFaceAlpha = 1.0;
% scatter(x_well2, y_well2, 200, 'r', 'filled', 'HandleVisibility', 'off', 'MarkerFaceAlpha', 1.0);
% Overlay wells using scatter with 80% transparency
hWell1 = scatter(x_well1, y_well1, 300, 'r', 'filled', 'DisplayName', 'Well');
hWell1.MarkerFaceAlpha = 0.2;
scatter(x_well2, y_well2, 300, 'r', 'filled', 'HandleVisibility', 'off', 'MarkerFaceAlpha', 0.2);

% % Overlay river as a patch with NO transparency
% hRiver = patch(river_x, river_y, 'b', 'FaceAlpha', 1.0, 'EdgeColor', 'none', 'DisplayName', 'River');

% Overlay river as a patch with 80% transparency (FaceAlpha = 0.2)
hRiver = patch(river_x, river_y, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'River');

% Adjust ticks and add legend
xlim([0, NCOL*dx]);
ylim([0, NROW*dx]);  % Note: y-axis limit uses NCOL*dx? Ensure proper dimension.
xticks(0:dx:(NCOL*dx));
yticks(0:dy:(NROW*dy));
% Set font properties for the axes
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
% Create and customize the legend
lgd1 = legend([hWell1, hRiver], 'Orientation', 'horizontal', 'Location', 'southoutside');
lgd1.FontSize = 12;
lgd1.FontName = 'Times New Roman';

% 4) Plot 2: Filled Contour Plot (contourf)

figure;
contourf(Xc, Yc, H, 20, 'LineColor', 'none');
colormap('summer');
cb = colorbar;
cb.Label.String = 'Groundwater Head (ft)';
cb.FontSize = 12;
cb.FontName = 'Times New Roman';
cb.Label.FontSize = 12;
cb.Label.FontName = 'Times New Roman';
xlabel('x-coordinate (ft)');
ylabel('y-coordinate (ft)');
%title('A. Upper aquifer under natural conditions');
%title('B. Lower aquifer under natural conditions');
% title('C. Upper aquifer under pumping conditions ');
 title('D. Lower Aquifer under pumping conditions ');
axis equal;
grid on;
hold on;

% % Overlay wells using scatter with NO transparency
% hWell1 = scatter(x_well1, y_well1, 200, 'r', 'filled', 'DisplayName', 'Well');
% hWell1.MarkerFaceAlpha = 1.0;
% scatter(x_well2, y_well2, 200, 'r', 'filled', 'HandleVisibility', 'off', 'MarkerFaceAlpha', 1.0);
% Overlay wells using scatter with 80% transparency
hWell1 = scatter(x_well1, y_well1, 300, 'r', 'filled', 'DisplayName', 'Well');
hWell1.MarkerFaceAlpha = 0.2;
scatter(x_well2, y_well2, 300, 'r', 'filled', 'HandleVisibility', 'off', 'MarkerFaceAlpha', 0.2);

% % Overlay river as a patch with NO transparency
% hRiver = patch(river_x, river_y, 'b', 'FaceAlpha', 1.0, 'EdgeColor', 'none', 'DisplayName', 'River');
% Overlay river as a patch with 80% transparency
hRiver = patch(river_x, river_y, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'River');

% Adjust axes limits and ticks based on cell centers
xlim([min(x_centers)-dx/2, max(x_centers)+dx/2]);
ylim([min(y_centers)-dy/2, max(y_centers)+dy/2]);
xticks(0:dx:(NCOL*dx));
yticks(0:dy:(NROW*dy));
% Set font properties for the axes
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

% Create and customize the legend
lgd2 = legend([hWell1, hRiver], 'Orientation', 'horizontal', 'Location', 'southoutside');
lgd2.FontSize = 12;
lgd2.FontName = 'Times New Roman';


