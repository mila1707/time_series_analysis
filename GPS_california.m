% The analysis of the PBO network
clc; clear all; close all;


% add the Code to the path
addpath(genpath('./'));
% list the stations in the folder
!ls -1 ../data_nam08_orig/*.csv > stations.list

%% read the data
fprintf('importing GPS data ... \n');
GPS_dataset = 'stations.list';
GPS_fid_all = fopen(GPS_dataset, 'r');
while ~feof(GPS_fid_all)
%     for ii=1:10
ii=1;

% compute the length of the rows
fid=fopen('stations.list','rt'); % t???fread???????
row=0;
while ~feof(fid)
    % ?????10000???????????????10????ASCII??
    % '*char'???????????*????????
    % ??fread???????????????????????
    % ??fopen????????????gbk
    row=row+sum(fread(fid,10000,'*char')==char(10));
    % ????????????????????????????
    % 'char'?????????????????double?
    % ??????char????double?????????
    % row=row+sum(fread(fid,10000,'char')==10);
end
fclose(fid);
n=row;

GPS_input = fgetl(GPS_fid_all);
    
%% load data for test
GPS_input = '../data_nam08_orig/AB01.pbo.nam08.csv';
[GPS_Y,GPS_M,GPS_D,GPS_N,GPS_E,GPS_U,GPS_Ns,GPS_Es,GPS_Us,GPS_Qu,GPS_DecYr] = load_data(GPS_input);

% GPS_input2 = '../data_nam08_orig/P270.pbo.detrended.nam08.csv';
% [GPS_Y2,GPS_M2,GPS_D2,GPS_N2,GPS_E2,GPS_U2,GPS_Ns2,GPS_Es2,GPS_Us2,GPS_Qu2,GPS_DecYr2] = load_data(GPS_input2);

%% interpolation of GPS time series
range1=2008;
range2=2018;
[GPS_Ni2,GPS_Ei2,GPS_Ui2,decdatei2] = GPS_interpolation(GPS_DecYr,GPS_N,GPS_E,GPS_U,range1,range2);

%% detrend the time series
GPS_dN=detrend_GPS(GPS_DecYr,GPS_N);
trend_N=GPS_N-GPS_dN;
GPS_dE=detrend_GPS(GPS_DecYr,GPS_E);
trend_E=GPS_E-GPS_dE;
GPS_dU=detrend_GPS(GPS_DecYr,GPS_U);
trend_U=GPS_U-GPS_dU;


%% save the time series
% GPS_DecYrM=zeros(length(GPS_DecYr),n);
% % GPS_siteM=table(length(GPS_DecYr),3)
% GPS_DecYrM(:,ii)=GPS_DecYr;
% % GPS_siteM(:,1)=GPS_site;

%% fit the time series
fun=@(x,xdata) x(1)+x(2)*xdata+x(3)*sin(2*pi*xdata+x(4)/180*pi)+x(5)*sin(pi*xdata+x(6)/180*pi)+x(7)*heaviside(xdata-0);
[GPS_U_fit,resnorm,residual,exitflag,output,lambda,jacobian] = ...
        lsqcurvefit(fun,...
        [0.0001,0.000001,0.00001,0.00001,0.00001,0.00001,0.00001],GPS_DecYr,GPS_dU);
[GPS_N_fit,resnorm,residual,exitflag,output,lambda,jacobian] = ...
        lsqcurvefit(fun,...
        [0.0001,0.000001,0.00001,0.00001,0.00001,0.00001,0.00001],GPS_DecYr,GPS_dN);
[GPS_E_fit,resnorm,residual,exitflag,output,lambda,jacobian] = ...
        lsqcurvefit(fun,...
        [0.0001,0.000001,0.00001,0.00001,0.00001,0.00001,0.00001],GPS_DecYr,GPS_dE);
    
    
%% Plot the time series
% !psxy 
figure;
subplot(3,1,1);
plot(GPS_DecYr,GPS_N,'.r');
hold on
plot(GPS_DecYr,trend_N,'.g','linewidth',0.5);
plot(decdatei2,GPS_Ni2,'.k');
yyaxis right;
plot(GPS_DecYr,GPS_dN,'.b');
% plot(GPS_DecYr,fun(GPS_N_fit,GPS_DecYr),'-k','linewidth',2);
% plot(GPS_DecYr2,GPS_N2,'.k');
hold off
ylabel('displacement(mm)');
% set the title of the figure
title(GPS_input(20:23));
set(gca,'XTick',round(min(GPS_DecYr)):1:round(max(GPS_DecYr)));
subplot(3,1,2);
plot(GPS_DecYr,GPS_E,'.r');
hold on
plot(GPS_DecYr,trend_E,'.g','linewidth',0.5);
plot(decdatei2,GPS_Ei2,'.k');
yyaxis right;
plot(GPS_DecYr,GPS_dE,'.b');
% plot(GPS_DecYr,fun(GPS_E_fit,GPS_DecYr),'-k','linewidth',2);
% plot(GPS_DecYr2,GPS_E2,'.k');
% plot(decdatei2,GPS_Ei2,'.k');
hold off
ylabel('displacement(mm)');
set(gca,'XTick',round(min(GPS_DecYr)):1:round(max(GPS_DecYr)));
subplot(3,1,3);
plot(GPS_DecYr,GPS_U,'.r');
hold on
plot(GPS_DecYr,trend_U,'.g','linewidth',0.5);
plot(decdatei2,GPS_Ui2,'.k');
yyaxis right;
plot(GPS_DecYr,GPS_dU,'.b');
% plot(GPS_DecYr,fun(GPS_U_fit,GPS_DecYr),'-k','linewidth',2);
% plot(GPS_DecYr2,GPS_U2,'.k');
hold off
xlabel('epoch(yr)');
ylabel('displacement(mm)');
set(gca,'XTick',round(min(GPS_DecYr)):1:round(max(GPS_DecYr)));

ii=ii+1;

break;
end





