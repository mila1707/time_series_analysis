function date_time = yrdoy2date(year,doy)

% yrdoy2date.m converts year and doy into date (year, month, day, hours,
% minutes, seconds).
% -------------------------------------------------------------------------
% INPUT
% year: Vector containing the year. Size: Tx1
% doy: Vector containing the day of the year. Size: Tx1
% -------------------------------------------------------------------------
% OUTPUT
% date_time: Matrix containing the year, month, day, hours, minutes, and
%            seconds. Size: Tx6
% -------------------------------------------------------------------------
% 
% Adriano Gualandi - 21 Mar 2017
% California Institute of Technology
% Geological and Planetary Science Division

% Find the indexes for which the year is a leap year
%ind_leap = leapyear(year);
ind_leap = leapyear_notlicensed(year);
% Create two lists of days for regular and leap years
days      = [31,28,31,30,31,30,31,31,30,31,30,31];
days_leap = [31,29,31,30,31,30,31,31,30,31,30,31];
% Create two lists of cumulative days for regular and leap years
cum_days      = cumsum([0,days]);
cum_days_leap = cumsum([0,days_leap]);
% Find the total number of dates to be estimated
T = numel(year);
% Reshape the input variable year if it has a size 1xT instead of Tx1
if size(year,2) == T
    year = year';
end
% Check that the input variables have the same size
if size(year)~=size(doy)
    % If the variables do not have the same size, transpose the second
    % input variable doy
    doy = doy';
end
% Check again that the input variables have the same size
if size(year)~=size(doy)
    % If not, then there is an error
    error('The input variables year and doy must be vectors of the same size.');
end
% For every date to be estimated, find the indexes where the count of
% cumulative days exceeds the doy. The ind_month_all variable is a logical
% variable that is equal to 0 until the desired month(+1) is reached.
ind_month_all = (repmat(doy,1,13)-repmat(cum_days,T,1))<0;
% Find the indexes for which the first 1 is encountered in the logical
% matrix ind_month_all
[~,ind_month] = max(ind_month_all,[],2);
% The month is the index minus 1 because the cum_days variable contains 13
% entries
month = ind_month-1;
% Repeat the same but for leap years
ind_month_leap_all = (repmat(doy,1,13)-repmat(cum_days_leap,T,1))<0;
[~,ind_month_leap] = max(ind_month_leap_all,[],2);
% Substitute the month at the ind_leap positions
month(ind_leap) = ind_month_leap(ind_leap)-1;
month(month==0) = 1;
% Calculate the day subtracting the cumulative number of days of the month
% calculated at previous step from the day of the year
day_decimal = doy-cum_days(month)';
day_decimal(ind_leap) = doy(ind_leap)-cum_days_leap(month(ind_leap))';
day = 1+fix(day_decimal);
% The remainder is the hours, minutes, seconds in decimal format
hrminsec = day_decimal-fix(day_decimal);
% Find the hours
hours = fix(hrminsec*24);
% The reminder is the minutes, seconds in decimal format
minsec = hrminsec*24 - hours;
% Find the minutes
minutes = fix(minsec*60);
% The reminder is the seconds in decimal format, that multiplied by 60
% gives the seconds
seconds = (minsec*60-minutes)*60;
% The output variable is the matrix containing the year, month, day, hours,
% minutes, and seconds variables
date_time = [year,month,day,hours,minutes,seconds];