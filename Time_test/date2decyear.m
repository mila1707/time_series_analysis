function timeline = date2decyear(date_time)
% date2decyear.m converts time from date format to decimal year format.
% -------------------------------------------------------------------------
% INPUT
% date_time: Matrix containing the year, month, day, hours, minutes, and
%            seconds. Size: Tx6
%            Format:
%             Column 1: year
%             Column 2: month
%             Column 3: hours
%             Column 5: minutes
%             Column 6: seconds
% -------------------------------------------------------------------------
% OUTPUT
% timeline: Vector containing the date in decimal year. Size: Tx1
% -------------------------------------------------------------------------
% 
% Jun Yan - 20 Mar 2019
% UC Berkeley
% Earth and Planetary Science Department

% Initialize the starting variables
year    = date_time(:,1);
month   = date_time(:,2);
day     = date_time(:,3);
hours   = date_time(:,4);
minutes = date_time(:,5);
seconds = date_time(:,6);

% Create two lists of cumulative days for regular and leap years
cum_days = [0,31,59,90,120,151,181,212,243,273,304,334]';
cum_days_leap = [0,31,60,91,121,152,182,213,244,274,305,335]';

% Calculate the number of the date
[num_date,~]=size(date_time);
% Convert the date format in decimal year format
% Check if the year is regular or leap
for i=1:num_date
% If the date belongs to a regular year...
    if ~leapyear_notlicensed(year(i))
% Convert the month and day in the day of the year (doy)
        doy(i,1)=cum_days(month(i))+day(i)-1;
% Calculate the decimal year by using the variables of year and doy
        timeline(i,1)=year(i,1)+doy(i,1)./365;
% If the date belongs to a leap year...
    else
        doy(i,1)=cum_days_leap(month(i))+day(i)-1;
        timeline(i,1)=year(i,1)+doy(i,1)./366;
    end
end
% Calculate the decimal year by using the variables of year and doy
% % If the date belongs to a regular year...
% if ~leapyear_notlicensed(year
% timeline=year+doy./yeardays(year);

end

