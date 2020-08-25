% File name: california_test.m
% Author: Jun Yan
% Version: 1.1
% Date: 12/25/2019
% modified: 7/13/2020
% Description: This code is for testing the Regam. Here we try to use the
% Regam to replace missing values in GPS time series with imputed values.  

clc;close all;
clear all;

% STEP 1
% make a list of the file name
!ls -1 ../data_nam08_test/*.csv > GPS.filename
% read the name of the file
filename=cell2mat(textread('GPS.filename','%s'));

% assume that X is a matrix with n*p, which n is the epoch of observation,
% and p is the number of variance
% initialize the input matrix X2

% set the range of date by the user first
begin_year = 2011.0;
begin_month = 1;
begin_day = 1;
end_year = 2017.0;
end_month = 1;
end_day =1;
span_years = [begin_year:1:end_year];

% read the GPS file and extract the data
[num_file,~]=size(filename);
for ii=1:num_file
    file_id=fopen(filename(ii,:));
    disp(['reading ',filename(ii,:)]);
    s=dir(filename(ii,:));
    if s.bytes == 0
        disp('this file is empty');
        continue;
    else
        X0=textscan(file_id,'%d-%d-%d %f %f %f %f %f','delimiter',',');
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change the string of date to decimal date
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        days_to_month=[0,31,59,90,120,151,181,212,243,273,304,334]';
        days_to_month_leapyear=[0,31,60,91,121,152,182,213,244,274,305,335]';
        [num_date,~]=size(X0{1});
        year1=X0{1};
        month1=X0{2};
        day1=X0{3};
        doy=[];
        decdate=[];
        for i=1:num_date
            doy(i,1)=days_to_month(month1(i))+day1(i)-1;
            if leapyear(year1(i))
                doy(i,1)=days_to_month_leapyear(month1(i))+day1(i)-1;
            end
        end
        decdate=double(year1)+double(doy)./yeardays(double(year1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % reprint the data by using the decdate format which including the date and
        % E N U deformation
        X1=[decdate X0{4} X0{5} X0{6}];

% find the position of missing data and fill with NaN
        ind = 1;
        ind_filled=1;
        for year=span_years(1):span_years(end)
            if leapyear(year)
                days = 366;
            else
                days = 365;
            end
            for i=1:days
                standard_decdate = double(year)+i*(1.0/days);
                epoch(ind_filled)=standard_decdate;
                while 1 && ind < length(decdate)
                    if abs(decdate(ind)-standard_decdate) < 0.0002
                        X_NaN(ind_filled,:)=X1(ind,:);
                        ind_filled = ind_filled +1;
                        ind = ind +1;
                        break;
                    elseif decdate(ind)-standard_decdate > 0.0025
                        X_NaN(ind_filled,:)=[standard_decdate NaN NaN NaN];
                        ind_filled = ind_filled + 1;
                        break;
                    elseif standard_decdate-decdate(ind) > 0.0025
                        ind = ind +1;
                    end
                end
            end
        end                

        % assemble the X_filled to X_input
        if ii==1
%             X_input_N=[X_NaN(:,3)];
%             X_input_V=[X_NaN(:,4)];
            X_input=[X_NaN(:,4)];
        else
            X_input=[X_input X_NaN(:,4)];
        end
    end
    fclose(file_id);
end

% RUN the "regem" to interpolate the data X
%     Field name         Parameter                                  Default
% OPTIONS.regress    Regression procedure to be used:           'mridge'
%                        'mridge': multiple ridge regression
%                        'iridge': individual ridge regressions
%                        'ttls':   truncated total least squares 
%                                  regression 
% DEFINE the varible OPTIONS
OPTIONS.regress='mridge';

addpath('/data/home/jyan/MATLAB/RegEM-master');
[X, M, C, Xerr, B, peff, kavlr, kmisr, iptrn] = regem(X_input, OPTIONS);

% save the interpolated data
for iii=1:num_file
%     X_filled=[epoch(:)',X(:,iii)];
    fileid=fopen([filename(iii,:),'_filled'],'wt');
    for epoch_i=1:length(epoch)
        fprintf(fileid,'%9.6f, ',epoch(epoch_i)');
        fprintf(fileid,'%9.2f\n',X(epoch_i,iii));
    end
end   
        
    
% figure;
% plot(X_NaN(:,1),X_NaN(:,4),'.-');
% figure;
% plot(epoch(:)',X(:,98),'.-');

% % make a figure of X
% figure;
% plot(X(:,1),X(:,4),'.-');

