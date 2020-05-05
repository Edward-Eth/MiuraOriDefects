close all

%If the matrix hasn't already been loaded in from the spreadsheets then do
%so, if it has then don't do it again because it takes ages. To perform the
%matrix read command.
if ~exist('Matrix','var') 
    clear all

    %Fetch all excel files in this directory
    Files=dir('*.xls');
    for k=1:length(Files)
       FileNames=Files(k).name;
    end

    %Save all fthe file names in "SpreadSheet" for later use
    for ii=1:length(Files)
        SpreadSheet{ii} = Files(ii).name;
    end

    %Find out what kind of sweep run this is
    if contains(SpreadSheet{1},'method2')
        SweptVariable = 'Vertex Deviation';
    elseif contains(SpreadSheet{1},'method3')
        SweptVariable = 'Tear Proportion';
    elseif contains(SpreadSheet{1},'method5')
        SweptVariable = 'Fold Stiffness Variation';
    end
        
    
    %Fetch all data into Matrix and edit file names to more presentable
    %form
    for ii=1:length(SpreadSheet)
        Matrix{ii} = readmatrix(SpreadSheet{ii});
        SpreadSheet{ii} = strrep(SpreadSheet{ii},'_',' ');
        SpreadSheet{ii} = strrep(SpreadSheet{ii},'.xls','');
        SpreadSheet{ii} = strrep(SpreadSheet{ii},'runfolder ','');
        SpreadSheet{ii} = strrep(SpreadSheet{ii},' rep8','');
        SpreadSheet{ii} = strrep(SpreadSheet{ii},'method2',SweptVariable);
        SpreadSheet{ii} = strrep(SpreadSheet{ii},'method3',SweptVariable);
        SpreadSheet{ii} = strrep(SpreadSheet{ii},'method5',SweptVariable);
        SpreadSheet{ii} = strrep(SpreadSheet{ii},'dev','');
        SpreadSheet{ii} = strrep(SpreadSheet{ii},'tears','');
        SpreadSheet{ii} = strrep(SpreadSheet{ii},'stiffvary','');
        if strcmp(SpreadSheet{ii},'Tear Proportion')
            SpreadSheet{ii} = 'Untorn Control Case';
        end
        if strcmp(SpreadSheet{ii},'Fold Stiffness Variation')
            SpreadSheet{ii} = 'Uniform Stiffness Control Case';
        end
        if strcmp(SpreadSheet{ii},'Vertex Deviation')
            SpreadSheet{ii} = 'Undeviated Control Case';
        end
        VariableVal{ii} = strrep(SpreadSheet{ii},'Tear Proportion ','');
        VariableVal{ii} = strrep(VariableVal{ii},'Fold Stiffness Variation ','');
        VariableVal{ii} = strrep(VariableVal{ii},'Vertex Deviation ','');
        VariableVal{ii} = strrep(VariableVal{ii},'Untorn Control Case','');
        VariableVal{ii} = strrep(VariableVal{ii},'Uniform Stiffness Control Case','');
        VariableVal{ii} = strrep(VariableVal{ii},'Undeviated Control Case','');
        if strcmp(VariableVal{ii},'')
            VariableVal{ii} = '0';
        end
        VariableVal{ii} = str2double(VariableVal{ii});
    end
end

%Matrix {1,n} is all of the data from the nth run, contains (1,m) where m
%is the number of repeats run.

clearvars -except 'SpreadSheet' 'Matrix' 'VariableVal' 'SweptVariable'
%Clears all of the variables defined below this point at the start of each
%run to make sure there's no cross fuckery going on.

%Across all of the matrix
for ii=1:length(Matrix)
    clear 'startcolumns'
    clear 'endcolumns'
    clear 'starts'
    clear 'a'
    a = Matrix{1,ii}; %read in current property value's data
    
    %Each repeat is on the end of previous run, so must find the points in
    %the data where each run starts
    starts = zeros(1,length(a));
    starts(1) = 1;
    for ll=2:length(a)
        %start of a run is where the new x value is less than the previous,
        %more robost than using 10 as not all runs start at the same
        %values.
        if a(1,ll)<a(1,ll-1)
            starts(ll)=1;
        end
    end
    
    %Make a new list of the positions of the located start points
    startcolumns = find(starts==1);
    
    %Each run ends at the position before the start of the next, apart from
    %the last run which ends at the end of the table, added in the next
    %row.
    endcolumns = (startcolumns(2:end)-1);
    endcolumns = [endcolumns,length(a)];
    
    %Split up the data into seperate chunks using the starts and ends
    %found, as long as they're not an error (almost zero energy at all
    %points) run.
    for jj=1:(length(startcolumns)-1)
        if a(2,startcolumns(jj)+4)>0.001
            DataSets{ii,jj} = a(:,startcolumns(jj):endcolumns(jj));
        end
    end
    
    %Add the last one in outside of the loop for ease to prevent fence post
    %erroring
    jj=length(startcolumns);
    if a(2,startcolumns(jj)+5)>0.001
        DataSets{ii,jj} = a(:,startcolumns(jj):end);
    end
    
    %The above code sometimes leaves zeros in error runs instead of "empty"
    %cells, this code replaces all undesired data with emtpy cells for
    %later ease of use.
    for jj=1:length(DataSets(ii,:))
        if length(DataSets{ii,jj}) < 2
            DataSets{ii,jj} = [];            
        end
    end
end

for ii=1:length(DataSets)

    %Clear all of the data from previous runs for sanity's sake. Shouldn't
    %matter but oh well.
    clear 'XVal';
    clear 'YVal';
    clear 'x';
    clear 'y';
    clear 'logx';
    clear 'logy';
    clear 'Variance';
    
    %Don't clear CurrentData, just set it empty so you can concatenate in
    %the following loop
    CurrentData = [];
    
    %Recollect all non-error repeats into one data set for averaging
    %purposes. This effectively returns the data to the intitial state
    %minus the errored runs.
    for jj=1:length(DataSets(ii,:))
        if ~isempty(DataSets(ii,jj))
            CurrentData = [CurrentData,DataSets{ii,jj}];
        end
    end
    
    %Set up limits on the figure so that the best fit curve knows how much
    %of x to plot over.
    figure
    hold on
    xlim([0 100])
    ylim([0 inf])
    
    %Sort all the data into ascending x order
    fittingData = sortrows(CurrentData');
    
    %Initial values
    currentX = fittingData(1,1);
    currentXcount = 0;
    currentYsum = 0;
    xValTally = 1;
    
    for oo=1:length(fittingData(:,1))
       currentY = fittingData(oo,2);
       
       %If the new x value is equal to the "current" running X value, and
       %the y value is a number: add one to the number of data points
       %at the current x value, and add the current y value to the running
       %y total at this x point.
       if fittingData(oo,1) == currentX && ~isnan(fittingData(oo,2))
           currentXcount = currentXcount+1;
           currentYsum = currentYsum + currentY;
           YValuesUsed(xValTally,currentXcount) = currentY;
           
       %If the new x value is not the same as the running X value, there
       %were more than three data points at the previous x value and the
       %current y value is a number:
       elseif currentXcount > 2 && ~isnan(fittingData(oo,2))
           YCalc = currentYsum/currentXcount;
           
           %If the variable YVal exists:
           if exist('YVal','var')
               
                %If the change between the new Yvalue to be added and the
                %previous is less than 0.02: set the new YValue to be the
                %average of the running y value total at the x point,
                %record how many data points were averaged to get this y
                %value and store the x value.
                if abs(YCalc-YVal(end))<0.01
                    YVal(xValTally)=YCalc;
                    Variance(xValTally) = var(nonzeros(YValuesUsed(xValTally,1:6)));
                    AveragePointCount(xValTally) = currentXcount;
                    XVal(xValTally)=currentX;
                    
                %If the change was more than 0.02, reset the count of the
                %number of x points so that there is not an empty column.
                else
                    xValTally = xValTally-1;
                end
                
           %If YVal doesn't exist alread, store the y,x and contributing point
           %count without checking them.
           else
               YVal(xValTally)=YCalc;
               Variance(xValTally) = var(nonzeros(YValuesUsed(xValTally,:)));
               AveragePointCount(xValTally) = currentXcount;
               XVal(xValTally)=currentX;
           end
           
           %Set running x to the next xvalue in the list of x values, add
           %one to the count of how many x values there have been, reset
           %the running total values
           currentX = fittingData(oo+1,1);
           xValTally = xValTally+1;
           currentXcount=0;
           currentYsum=0;
           
       %Else if the current oo count is less than the length of list of
       %values: set the current x to be the next entry. (This prevents
       %overrun of the values matrix.
       elseif oo < length(fittingData(:,1))
           currentX = fittingData(oo+1,1);
           currentXcount=1;
           currentYsum=currentY;
       end
    end
    
    %Rename variables for ease.
    x=XVal';
    y=YVal';
    
    %Extract the variance at x=50 for comparison purposes
    loc=find(x==50);
    Variance50 = Variance(loc);
    
    %Take the log of the x and y values to find the correct power to raise
    %x to for the best fit curve.
    logx = log(x);
    logy = log(y);
    
    %Linear fit to log values to find n in x^n for fitting the actual data.
    [PowerFinder,~] = fit(logx,logy,'poly1');
    Coeffs = coeffvalues(PowerFinder);
    Power = Coeffs(1);
    
    %Create custom fit string with the power value that was just found.
    fitstring = strcat('x^',num2str(Power));
    testfit = fittype({fitstring});
    
    %Fit the current average data using the above fit.
    [fitcurve{ii},gof{ii}] = fit(x,y,testfit);
    
    %Export the important data found to one matrix for plotting and
    %comparison
    Multiplier = coeffvalues(fitcurve{ii});
    Parameters{ii,1} = Power;
    Parameters{ii,2} = Multiplier;
    Parameters{ii,3} = Power*Multiplier;
    Parameters{ii,4} = VariableVal{ii};
    Parameters{ii,5} = SpreadSheet{ii};
    Parameters{ii,6} = x;
    Parameters{ii,7} = y;
    Parameters{ii,8} = YValuesUsed;
    Parameters{ii,9} = Variance;
    Parameters{ii,10} = Variance50;
    
    %Set the graph title, axis labels etc
    title(SpreadSheet{ii})
    ylabel('Global Energy Absorbtion (J)')
    xlabel('Proportion of displacement')
    
    %Plot the line of best fit
    fitplot = plot(fitcurve{ii});
    set(fitplot,'LineWidth',3);
    set(fitplot,'DisplayName','Power Curve Trendline');
    
    %Place the legend
    legend('Location','best')
    
    %Set the limits to the actual desired ranges now that the best fit
    %curve has been plotted from x=0 to 100.
    xlim([0 105])
    ylim([0 inf])
    
    %Plot each valid data set on the same axis and set their legend entry
    %names as appropriate
    for kk=1:length(DataSets(ii,:))
        CurrentData=DataSets{ii,kk};
        if ~isempty(CurrentData) && length(CurrentData) ~= 1
            plot(CurrentData(1,:),CurrentData(2,:),'DisplayName',['Repeat',' ',num2str(kk)])
        end
    end
    
    %Plot the averaged data used for graph fitting
    plot(x,y,'DisplayName','Averaged Data')
    
    %Place the legend in top left corner
    legend('show','Location','northwest')
    
    %Save this figure to a file
    %saveas(gcf,strrep([strrep(SpreadSheet{ii},'.','_'),' ','Stiffness'],' ','-'),'pdf')
    print(strrep([strrep(SpreadSheet{ii},'.','_'),' ','Stiffness'],' ','-'),'-dpng','-r600')
    
    %Make a new figure and plot the fit residuals on them to check how good
    %the fit model is.
    figure
    plot(fitcurve{ii},x, y,'residuals');
    title([SpreadSheet{ii},' ','Fit Residuals'])
    %saveas(gcf,strrep([strrep(SpreadSheet{ii},'.','_'),' ','Fit Residuals'],' ','-'),'pdf')
    print(strrep([strrep(SpreadSheet{ii},'.','_'),' ','Fit Residuals'],' ','-'),'-dpng','-r600')
    
end

Parameters = sortrows(Parameters,4);
ValueNames = {'Power','Multiplier','Steepness Coefficient','Swept Variable Value','Plot Title','X Points for average curve','Y points for average curve','Y Values used to find average','Variance','Variance50'};
ValueTable = table(Parameters(:,1),Parameters(:,2),Parameters(:,3),Parameters(:,4),Parameters(:,5),Parameters(:,6),Parameters(:,7),Parameters(:,8),Parameters(:,9),Parameters(:,10),'VariableNames',ValueNames);

figure
title(['Steepness Parameter vs',' ',SweptVariable])
hold on
x = cell2mat(Parameters(:,4));
y = cell2mat(Parameters(:,3));
xlim([0 (1.1*x(end))])
ylim([0 inf])
Trendfit = fit(x,y,'poly2');
trendcurveplot = plot(Trendfit);
set(trendcurveplot,'DisplayName','Quadratic Trend Line');
ylabel('Stiffness Parameter')
xlabel(SweptVariable)
plot(x,y,'LineWidth',3,'DisplayName','Stiffness')
legend('Location','northwest')
%saveas(gcf,strrep(['Steepness Parameter vs',' ',SweptVariable],' ','-'),'pdf')
print(strrep(['Steepness Parameter vs',' ',SweptVariable],' ','-'),'-dpng','-r600')

% figure
% title('Variance of sample stiffnesses against displacement')
% hold on
% for ii=1:2:length(Parameters(:,1))
%     plot(cell2mat(ValueTable{ii,'X Points for average curve'}),cell2mat(ValueTable{ii,'Variance'}),'DisplayName',ValueTable{ii,'Plot Title'}{1})
% end
% xlim([0 100])
% ylim([0 inf])
% xlabel('Proportion of displacement')
% ylabel('Variance')
% legend('Location','northwest')

figure
title(['Stiffness Variance vs',' ',SweptVariable])
hold on
plot(cell2mat(ValueTable{1:20,'Swept Variable Value'}),cell2mat(ValueTable{1:20,'Variance50'}))
ylabel('Variance of Stiffness Values')
xlabel(SweptVariable)
%saveas(gcf,strrep(['Stiffness Variance vs',' ',SweptVariable],' ','-'),'pdf')
print(strrep(['Stiffness Variance vs',' ',SweptVariable],' ','-'),'-dpng','-r600')