%close all
close all
clear
%load all variables
filepath = ('/Users/mkindel/Documents/MSanalysis_temp/2023_04_25_mk234_ss3_HIIT_cut_mat')
addpath(filepath)
%% 

load ('C.mat');
T = data; % time traces, unit by frame

%% 

load ('proj.mat');
proj = data;
load('A.mat');
I = data; % unit images, unit by row pixel by column pixel
%save("spatial_footprints_2.mat", "I", '-mat')
labels = coords.unit_labels.data';
bigI = I.*100000; %multiple each value to get a bigger range for image
% [h, l, w] = size(I);
%% 

F = 1:length(T); %set variable 1 to last frame
b = 8; % number of images in each round 
%rounds = ceil(h/4); % number of rounds, rounded up to the nearest integer
prompt = 'which units do you want to eliminate? e.g. [1 3 4] ';
keep =[];
load ('C.mat');
raw_trace = data; %for raw
%% 

clear data attrs coords dims name

[m,p] = size(T);
[k,l,o] = size(I);

% indices = find(labels(:,1) == -1);
% bestT = T;
% bestT(indices,:)=[];
% bestI = bigI;
% bestI(indices,:,:)=[];
% mergedI = merge(bestI);
% image(mergedI)

new_labels = zeros(length(labels), 1);
group = -1;
group_index = 1;
while group <= max(labels)
    new_labels(labels==group) = group_index;
    if group == max(labels)
        break
    else
        group = min(labels(labels>group));
        group_index = group_index + 1;
    end
end

num_groups = max(new_labels);

merged_units = ones (num_groups,p);
merged_units_I = ones(num_groups,l,o);

for i = 1:num_groups
    merged_units (i,:) = mean(T(new_labels == i,:),1);
    merged_units_I(i,:,:) = mean(bigI(new_labels == i,:,:),1);
end

merged_units_final = merged_units((2:end),:);
merged_units_final_I = merged_units_I((2:end),:,:);

for i = 1:num_groups
    merged_units_raw (i,:) = mean(raw_trace(new_labels == i,:),1);
end

merged_units_raw_final = merged_units_raw((2:end),:);

[m,n] = size(merged_units_final);

%% z-score fitted

meanVal = mean(merged_units_final,2);
stdVal = std(merged_units_final,[],2);

% meanVal = mean(merged_units_final(:,1:1000),2);
% stdVal = std(merged_units_final(:,1:1000),[],2);
% 
for i = 1:m
    for k = 1:n
        TminusBase(i,k) = merged_units_final(i,k) - meanVal(i);
        Tz(i,k) = TminusBase(i,k)./stdVal(i);
    end
end

%% z-score raw

meanValRaw = mean(merged_units_raw_final,2);
stdValRaw = std(merged_units_raw_final,[],2);

% meanVal = mean(merged_units_final(:,1:1000),2);
% stdVal = std(merged_units_final(:,1:1000),[],2);
% 
for i = 1:m
    for k = 1:n
        TminusBaseRaw(i,k) = merged_units_raw_final(i,k) - meanValRaw(i);
        TzRaw(i,k) = TminusBaseRaw(i,k)./stdValRaw(i);
    end
end


%% dF/F I think it works?

medianVal = median(merged_units_final(:,1:1000),2);

for i = 1:m
    for k = 1:n
        FminusMedian(i,k) = merged_units_final(i,k) - medianVal(i);
        dFF(i,k) = FminusMedian(i,k)./medianVal(i);
    end
end

%%
%mergedI = merge(merged_units_final_I);

% 
% for n = 1:length(labels)
%     if labels(n,1) ~= -1
%         for a = 1:
%         bestI(a,:,:) = I(n,:,:);
%         bestT(a,:) = T(n,:);
%         end
%     end
% end

% for k = 1:length(bestT)
%     if sum(bestT(k,2)) == 0
%         bestT(k) = [];
%     end
%     if sum(bestI(k,'all')) == 0
%         bestI(k) = [];
%     end
% end
%     
        


% for k = 1:rounds
%     
%     for i = 1:b
%         figure (1) % traces
%         subplot(4,2,i)
%         plot(F,T(i+(b*(k-1)),:))
%         title(['unit' num2str(i+(b*(k-1)))])
%     
%         figure (2) % images
%         subplot(4,2,i)
%         image(squeeze(bigI(i+(b*(k-1)),:,:)))
%         title(['unit' num2str(i+(b*(k-1)))])
%     
%     end
%     current_keep = input(prompt);
%     keep = [keep, current_keep];
% end

% newT = T(keep,:); %traces of all the kept units
% newI = bigI(keep,:,:); %matrix with images of all the kept units
% [m n] = size(TzRaw);
% sizeAdd = 5;
% for i = 1:m
%     newTplot(i,:) = Tz(i,:)+(i*sizeAdd);
% end
% 
% plot(F, newTplot,'LineWidth',1, Color="black")
% 
% 
% print('-painters','-depsc','myVectorFile')
