function [GPS_Y,GPS_M,GPS_D,GPS_N,GPS_E,GPS_U,GPS_Ns,GPS_Es,GPS_Us,GPS_Qu,GPS_DecYr] = load_data(GPS_input)
% load the data set of GPS time series
%   To get the ENU components and the decimal format of Date
%     GPS_input = fgetl(GPS_fid_all);
%     
% % GPS_input = '1LSU.pbo.nam08.csv';
% GPS_input = '../data_nam08_orig/P270.pbo.nam08.csv';
% GPS_input = 'AB01.pbo.nam08.csv';
GPS_fid   = fopen(GPS_input, 'r');
% read the data above '#'
count = 0;
while ~feof(GPS_fid)
    tline=fgetl(GPS_fid);
    if tline(1) ~= 'D'
        count = count+1;
    else
        break;
    end
end
GPS_C=textscan(GPS_fid,'%f %f %f, %f, %f, %f, %f, %f, %f, %s');
GPS_site = GPS_input(1:4);
GPS_Y = GPS_C{1};
GPS_M = -GPS_C{2};
GPS_D = -GPS_C{3};
GPS_N = GPS_C{4};
GPS_E = GPS_C{5};
GPS_U = GPS_C{6};
GPS_Ns = GPS_C{7};
GPS_Es = GPS_C{8};
GPS_Us = GPS_C{9};
GPS_Qu = GPS_C{10};

% change the string of date to decimal date
days_to_month=[0,31,59,90,120,151,181,212,243,273,304,334]';
days_to_month_leapyear=[0,31,60,91,121,152,182,213,244,274,305,335]';
[num_date,~]=size(GPS_Y);
year1=GPS_Y;
month1=GPS_M;
day1=GPS_D;
for i=1:num_date
    doy(i,1)=days_to_month(month1(i))+day1(i)-1;
    if leapyear(year1(i))
        doy(i,1)=days_to_month_leapyear(month1(i))+day1(i)-1;
    end
end
GPS_DecYr=year1+doy./yeardays(year1);

% % save the input data
% fid=fopen(save.txt
% clear doy tline
end

