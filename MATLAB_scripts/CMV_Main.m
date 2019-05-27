%% Load CSV files
file = 'C:/insertfilerdirectoryhere/synthetic_data2.csv';
data_raw = CMV_dataRead_v3(file);
data = data_raw{2};

% adjData = calibrateSensors(data);
% data(data<1) = 0;
% data = movingAverage(data,7);
len = length(data);
adjData = zeros(len,9);
for i = 1:9
    adjData_temp = sgolayfilt(data(:,i),2,3);
    adjData(:,i) = adjData_temp;
end
data = adjData;

% Normalize data
[row,col] = size(data);
data_norm = zeros(row,col);
for i = 1:col
    data_norm(:,i) = data(:,i)/max(abs(data(:,i)));
end

%% Get data sets
x_start = 318;%406500;
x_end = 331;%410000;

data_length = x_end - x_start + 1;
actual_sample = zeros(data_length,1,9);
actual_sample_norm = zeros(data_length,1,9);
for i = 1:9
    actual_temp = data(x_start:x_end,i);
    actual_norm_temp = data_norm(x_start:x_end,i);
    actual_sample(:,:,i) = actual_temp;
    actual_sample_norm(:,:,i) = actual_norm_temp;
end
actual_sample = reshape(actual_sample,data_length,9);
actual_sample_norm = reshape(actual_sample_norm,data_length,9);

%% Get matrices
cmv_sample = actual_sample; %specify data set to be analyzed
cmv_sample_norm = actual_sample_norm;
window = 10;
pages = data_length - window + 1;%ceil(data_length/window)+1

% Get raw pixels
I_raw = zeros(3,3,pages);
for j = 1:data_length
    I_temp = [cmv_sample(j,5) cmv_sample(j,4) cmv_sample(j,6);
              cmv_sample(j,8) cmv_sample(j,7) cmv_sample(j,9);
              cmv_sample(j,2) cmv_sample(j,1) cmv_sample(j,3)];
    
    I_raw(:,:,j) = I_temp;
end

% Get pixels
I = zeros(3,3,pages);
for j = 1:pages
    for k = 1:window
    I_temp = [cmv_sample(j-1+k,5) cmv_sample(j-1+k,4) cmv_sample(j-1+k,6);
              cmv_sample(j-1+k,8) cmv_sample(j-1+k,7) cmv_sample(j-1+k,9);
              cmv_sample(j-1+k,2) cmv_sample(j-1+k,1) cmv_sample(j-1+k,3)];
    
    I(:,:,j) = I(:,:,j) + I_temp;
    end
    I(:,:,j) = I(:,:,j)/window;
end

% Get normalized pixels
I_norm = zeros(3,3,pages);
for j = 1:pages
    for k = 1:window
    I_temp = [cmv_sample_norm(j-1+k,5) cmv_sample_norm(j-1+k,4) cmv_sample_norm(j-1+k,6);
              cmv_sample_norm(j-1+k,8) cmv_sample_norm(j-1+k,7) cmv_sample_norm(j-1+k,9);
              cmv_sample_norm(j-1+k,2) cmv_sample_norm(j-1+k,1) cmv_sample_norm(j-1+k,3)];
    
    I_norm(:,:,j) = I_norm(:,:,j) + I_temp;
    end
    I_norm(:,:,j) = I_norm(:,:,j)/window;
end

I_ave_norm = zeros(3,3,pages);
idx = 1;
for j = 1:data_length
    I_temp = [cmv_sample_norm(j,5) cmv_sample_norm(j,4) cmv_sample_norm(j,6);
              cmv_sample_norm(j,8) cmv_sample_norm(j,7) cmv_sample_norm(j,9);
              cmv_sample_norm(j,2) cmv_sample_norm(j,1) cmv_sample_norm(j,3)]
    
        I_ave_norm(:,:,idx) = I_ave_norm(:,:,idx) + I_temp;
    if mod(j,window) == 0
        I_ave_norm(:,:,idx) = I_ave_norm(:,:,idx)/window;
        idx = idx + 1;
    end
end

%% Find angle using Matrix Gradient Method
imData = I_ave_norm;
[image_row, image_col, ~] = size(imData);
theta = zeros(image_row,image_col,pages);
magnitude = zeros(image_row,image_col,pages);
theta2 = zeros(pages,1);
for idx = 1:(pages-1)
    [Gx, Gy] = imgradientxy(imData(:,:,idx),'sobel');    
    [Gmag, Gdir] = imgradient(Gx, Gy);
    
    theta(:,:,idx) = Gdir;
    magnitude(:,:,idx) = Gmag;     
    figure;quiver(Gx,Gy)
    MG(idx) = getframe(gcf);     
    close
end

% Find the tau for every minimum point
for i=1:9
[~,tau(i,1)]=min(cmv_sample_norm(:,i));
end

%% Plot Figures
figure
plot1 = plot(cmv_sample);
set(plot1(1),'DisplayName','South','LineWidth',2,'Color','k');
set(plot1(2),'DisplayName','South-West','LineWidth',2);
set(plot1(3),'DisplayName','South-East','LineWidth',2);
set(plot1(4),'DisplayName','Norfth','LineWidth',2);
set(plot1(5),'DisplayName','North-West','LineWidth',2);
set(plot1(6),'DisplayName','North-East','LineWidth',2);
set(plot1(7),'DisplayName','Origin','LineWidth',2);
set(plot1(8),'DisplayName','West','LineWidth',2);
set(plot1(9),'DisplayName','East','LineWidth',2);

%% Peak Matching
% Set which sensors to find the peaks and dips
sensor1 = cmv_sample_norm(:,2);
sensor2 = cmv_sample_norm(:,5);
sensor3 = cmv_sample_norm(:,3);
sensor4 = cmv_sample_norm(:,8);
sensor5 = cmv_sample_norm(:,1);

sensor1_inv=1./sensor1;
sensor2_inv=1./sensor2;
sensor3_inv=1./sensor3;
sensor4_inv=1./sensor4;
sensor5_inv=1./sensor5;

% Find the peaks
peak_distance = 60;
peak_prominence = 0.036;
peak_width = 15;

%% Plot peaks and dips
figure

[p1,k1] = findpeaks(sensor1,'MinPeakProminence',peak_prominence,'MinPeakDistance',peak_distance,'MinPeakWidth',peak_width);
[dp1, dk1] = findpeaks(sensor1_inv,'MinPeakProminence',peak_prominence,'MinPeakDistance',peak_distance,'MinPeakWidth',peak_width);
s1 = plot(t,sensor1,'DisplayName','Origin Sensor','LineWidth',2);
hold on;

plot(t(k1),p1,'ro','markersize',10);
plot(t(dk1),1./dp1,'rs','markersize',10);

[p2,k2] = findpeaks(sensor2,'MinPeakProminence',peak_prominence,'MinPeakDistance',peak_distance,'MinPeakWidth',peak_width);
[dp2, dk2] = findpeaks(sensor2_inv,'MinPeakProminence',peak_prominence,'MinPeakDistance',peak_distance,'MinPeakWidth',peak_width);
s2 = plot(t,sensor2,'DisplayName','North Sensor','LineWidth',2);

plot(t(k2),p2,'ro','markersize',10);
plot(t(dk2),1./dp2,'rs','markersize',10);

[p3,k3] = findpeaks(sensor3,'MinPeakProminence',peak_prominence,'MinPeakDistance',peak_distance,'MinPeakWidth',peak_width);
[dp3, dk3] = findpeaks(sensor3_inv,'MinPeakProminence',peak_prominence,'MinPeakDistance',peak_distance,'MinPeakWidth',peak_width);
s3 = plot(t,sensor3,'DisplayName','East Sensor','Color','m','LineWidth',2);

plot(t(k3),p3,'ro','markersize',10);
plot(t(dk3),1./dp3,'rs','markersize',10);

hold off;
legend([s1 s2 s3], 'Origin Sensor','North Sensor','East Sensor')

% circshift lookup of local minima

%% Find the optical flow
F1 = getOpticalFlow(imData);
