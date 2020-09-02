function MarioBioAnalysisv1(knn,numyears,spec) 
% I have this written as a function but if you'd just like to take and 
% use bits and pieces you can just delete the first line and the "end"
% at the bottom. If you'd like to just run the function but knn=0 if you'd
% like to use linear interpolation to fit missing data or knn > 0 to use
% k-nearest neighbors. Set numyears to be the number of years you'd like to
% inlude if doing linear interpolation and spec = 'P', 'S', or 'Both'
% depending on what species you'd like to look at

knn =0;
numyears = 30;
spec = 'Both';

%% Load Data 
numData = importdata('LinExt_Eric.csv'); %Loads the numerical data
metaData = array2table(table2array(readtable('LinExt_metaData.csv','ReadRowNames',true)).'); %Loads the metadata

%Separate the data from each of the cores from the years
cores = numData.data(:,2:end);
years = numData.data(:,1);

%Calculate the number of rows and columns in cores
[numRows numCores] = size(cores);

%% Show Core Lengths

%Create boolean vector of all missing values 
% 1 = not missing, 0 = missing
missing = ~isnan(cores); 

%Initialize lastObs vector
lastObs = zeros(numCores,1);

%Loop through each core and calculate the last observation which is not
%missing for each core
for i = 1:numCores
    lastObs(i) = find(missing(:,i),1,'last');
end

%sort the cores/metadata in decreasing order by how far back they go
[yearsSort indices] = sort(lastObs, 'descend');
coreSort = cores(:,indices);
metaSort = metaData(indices,:);

%The following code is used to generate the plot of core length
%Against the number of available cores
[empMzr, cutoffs] = ecdf(years(lastObs));
plot(cutoffs,empMzr*numCores,'LineWidth',2)
title('Core Lengths vs. Time','FontSize', 24)
xlabel('Year Available','FontSize', 18)
ylabel('Number of Cores','FontSize', 18)
hold on
plot([1984 1984],[0 70],'LineWidth',2,'Color',[1,0,0])
plot([1974 1974],[0 70],'LineWidth',2,'Color',[1,0,0])
plot([1964 1964],[0 70],'LineWidth',2,'Color',[1,0,0])
plot([1954 1954],[0 70],'LineWidth',2,'Color',[1,0,0])
fig = gcf;
set(gcf, 'Position', get(0, 'Screensize'));
print(fig,'CoreLengths','-dpng');
hold off


%% Handle Missing Data

%Then next bit of code basically computes how many year to include in the 
%analysis, how many cores to include in the analysis, and what method to
%use for filling in the missing values
% In general: obs = number of cores to use, numyears = numbers of years to
% use, and knn = how many nearest neighbors to use if using k-nearest
% neighbors
if knn == 0 %Use linear interpolation
    obs = sum(years(lastObs)<=years(numyears)); %Sets obs equal to the number of cores that go past numyears
    %Note that you will need the most recent version of MatLab to use the
    %next function
    coresFilledTemp = fillmissing(coreSort(1:numyears,1:obs),'linear'); %Linear interpolation
else %Use k-nearest neighbors
    numyears = 137; %include all years and observations
    obs = 67;
    coresFilledTemp = knnimpute(coreSort,knn); %Fill missing data
end

%The next bit of code separates out the species
if strcmp(spec,'P') || strcmp(spec,'S')
    metaTemp = metaSort(1:obs,:); %Include only observations that are long enough
    meta = metaTemp(ismember(metaTemp.Var1, spec),:); %include only metadata corresponding to correct species
    coresFilled = coresFilledTemp(:,ismember(metaTemp.Var1, spec)); %include core data corresponding to correst species
else
    meta = metaSort(1:obs,:); %Include only observations that are long enough
    coresFilled = coresFilledTemp; %Including both species
    spec = 'Both'; %This is used in file output
end

%Calculate the number of observations included. Important if you select a
%species
[~, obs] = size(coresFilled);

%The next chunk of code is used a lot of Marron's functions. Many times you
%need to put in seperate colors for any classes that you want to visually
%separate. Colors are coded with three numbers as [a b c] so thats whats
%going on
species = zeros(obs,3);
species(ismember(meta.Var1, 'P'),1) = 1;
species(ismember(meta.Var1, 'S'),3) = 1;
site = zeros(obs,3);
site(ismember(meta.Var2, 'CR'),1) = 1;
site(ismember(meta.Var2, 'BH'),2) = 1;
site(ismember(meta.Var2, 'WW'),3) = 1;
site(ismember(meta.Var2, 'BS'),[1 2]) = 1;
site(ismember(meta.Var2, 'AR'),[1 3]) = 1;
site(ismember(meta.Var2, 'ES'),[2 3]) = 1;
site(ismember(meta.Var2, 'CF'), 1) = .5;
site(ismember(meta.Var2, 'FR'), 2) = .5;
transect = zeros(obs,3);
transect(ismember(meta.Var3, 'UK1'),[1 3]) = 1;
transect(ismember(meta.Var3, 'UK2'),1) = 1;
transect(ismember(meta.Var3, 'UK3'),2) = 1;
transect(ismember(meta.Var3, 'LK'),3) = 1;
reefzone = zeros(obs,3);
reefzone(ismember(meta.Var4, 'IR'),1) = 1;
reefzone(ismember(meta.Var4, 'OR'),3) = 1;

reefzone2 = zeros(obs,1);
reefzone2(ismember(meta.Var4, 'IR')) = 1;
reefzone2(ismember(meta.Var4, 'OR')) = 2;

irs = coresFilled(:,reefzone2 ==1)
ors = coresFilled(:,reefzone2 ==2)
grpstats(coresFilled(17,:),reefzone2)


%The following chunk of code will the first few principal components for
%reefzone. Not the same as a PCA plot. This will display the actual PCA
%directions rather than a scatter plot of the data projected onto those
%directions (that's in the next section
paramstruct = struct('npcadiradd', 4, ... %Number of PCA directions you want
                        'icolor', reefzone(1:obs,:), ... %You can reperpose this but putting a different variable here
                        'isubpopkde', 1, ...
                        'legendcellstr', {{'IR','OR'}}, ... %Lables the classes
                        'mlegendcolor', [1 0 0; 0 0 1]); %Colors that you want
curvdatSM(coresFilled,paramstruct)
fig = gcf;
    set(gcf, 'Position', get(0, 'Screensize'));
    print(fig,strcat('RZPCA_Spec_',spec,num2str(numyears),'knn',num2str(knn)),'-dpng');
%% Last 30 Years PCA Analysis

%The next four plots are PCA plots for each of the four metadata variables
%and the parameters are similar in nature to the last function described
%above
paramstruct = struct('npcadiradd', 4, ...
                        'icolor', reefzone(1:obs,:), ...
                        'isubpopkde', 1, ...
                        'legendcellstr', {{'IR','OR'}}, ...
                        'mlegendcolor', [1 0 0; 0 0 1]);
scatplotSM(coresFilled,[],paramstruct)
fig = gcf;
    set(gcf, 'Position', get(0, 'Screensize'));
    print(fig,strcat('ReefZonePCA_Spec_',spec,num2str(numyears),'knn',num2str(knn)),'-dpng');

paramstruct = struct('npcadiradd', 4, ...
                        'icolor', site(1:obs,:), ...
                        'isubpopkde', 1, ...
                        'legendcellstr', {{'CR','BH','WW','BS','AR','ES','CF','FR'}}, ...
                        'mlegendcolor', [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; .5 0 0; 0 .5 0]);
scatplotSM(coresFilled,[],paramstruct)
fig = gcf;
    set(gcf, 'Position', get(0, 'Screensize'));
    print(fig,strcat('SitePCA_Spec_',spec,num2str(numyears),'knn',num2str(knn)),'-dpng');

paramstruct = struct('npcadiradd', 4, ...
                        'icolor', species(1:obs,:), ...
                        'isubpopkde', 1, ...
                        'legendcellstr', {{'P','S'}}, ...
                        'mlegendcolor', [1 0 0; 0 0 1]);
scatplotSM(coresFilled,[],paramstruct)
fig = gcf;
    set(gcf, 'Position', get(0, 'Screensize'));
    print(fig,strcat('SpeciesPCA_Spec_',spec,num2str(numyears),'knn',num2str(knn)),'-dpng');

paramstruct = struct('npcadiradd', 4, ...
                        'icolor', transect(1:obs,:), ...
                        'isubpopkde', 1, ...
                        'legendcellstr', {{'UK1','UK2','UK3','LK'}}, ...
                        'mlegendcolor', [1 0 1; 1 0 0; 0 1 0; 0 0 1]);
scatplotSM(coresFilled,[],paramstruct)
fig = gcf;
    set(gcf, 'Position', get(0, 'Screensize'));
    print(fig,strcat('TransectPCA_Spec_',spec,num2str(numyears),'knn',num2str(knn)),'-dpng');

%% Last 30 Years DWD Analysis

%Separate IR and OR cores into separate data sets
ir = coresFilled(:,ismember(meta.Var4, 'IR'));
or = coresFilled(:,ismember(meta.Var4, 'OR'));

%Calculate DWD Directions between IR and OR. You can certainly use the
%function for other classes if you'd like. Just replace 'ir' and 'or' with
%data from the two classes you're interested in
dwddir = DWD2XQ(ir,or);

%Similar to the PCA plots in the previous section but with the DWD
%Direction as the first direction
paramstruct = struct('npcadiradd', 3, ...
                        'icolor', reefzone(1:obs,:), ...
                        'isubpopkde', 1, ...
                        'legendcellstr', {{'IR','OR'}}, ...
                        'mlegendcolor', [1 0 0; 0 0 1]);
scatplotSM(coresFilled,dwddir,paramstruct)
fig = gcf;
    set(gcf, 'Position', get(0, 'Screensize'));
    print(fig,strcat('ReefzoneDWD_Spec_',spec,num2str(numyears),'knn',num2str(knn)),'-dpng');
 
%  %This function runs the DiProPerimTest
%  paramstruct = struct('npcadiradd', 3, ...
%                         'icolor',[1 0 0; 0 0 1], ...
%                         'isubpopkde', 1, ...
%                         'legendcellstr', {{'IR','OR'}}, ...
%                         'mlegendcolor', [1 0 0; 0 0 1]);
% [stat, pvalue, zscore] = DiProPermSM(ir,or,paramstruct);
% fig = gcf;
%     set(gcf, 'Position', get(0, 'Screensize'));
%     print(fig,strcat('ReefzoneDiProDerm_Spec_',spec,num2str(numyears),'knn',num2str(knn)),'-dpng');
% 
% clf
%Plot the loading on the DWD Direction
plot(years(1:numyears),dwddir)
fig = gcf;
    set(gcf, 'Position', get(0, 'Screensize'));
    print(fig,strcat('DWDDir_Spec_',spec,num2str(numyears),'knn',num2str(knn)),'-dpng');
    
clf  
plot(years(1:numyears),abs(dwddir))
fig = gcf;
    set(gcf, 'Position', get(0, 'Screensize'));
    print(fig,strcat('DWDDirAbs_Spec_',spec,num2str(numyears),'knn',num2str(knn)),'-dpng');
    
end