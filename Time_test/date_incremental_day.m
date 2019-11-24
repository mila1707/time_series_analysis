function date = date_incremental_day(start_date,daynum)

% date_incremental_day.m creates a list of dates given a starting date and
% the number of days after the starting date that the user desires.
% -------------------------------------------------------------------------
% INPUT
% start_date: Vector of size 1x6, containing the starting date. Format:
%             Column 1: year
%             Column 2: month
%             Column 3: day
%             Column 4: hours
%             Column 5: minutes
%             Column 6: seconds
% daynum: Number of days that you want to convert into date format after
%         the start_date. Scalar variable.
% -------------------------------------------------------------------------
% OUTPUT
% date: Matrix of size daynumx6, containing the list of the dates. Format:
%       Column 1: year
%       Column 2: month
%       Column 3: day
%       Column 4: hours
%       Column 5: minutes
%       Column 6: seconds
% -------------------------------------------------------------------------
% 
% Adriano Gualandi - 10 Feb 2017
% California Institute of Technology
% Geological and Planetary Science Division
% modified by Jun Yan - 21 Mar 2019

% Initialize the starting variables
year    = start_date(1);
month   = start_date(2);
day     = start_date(3);
hours   = start_date(4);
minutes = start_date(5);
seconds = start_date(6);
% Create two lists of days for regular and leap years
days = [31,28,31,30,31,30,31,31,30,31,30,31];
days_leap = [31,29,31,30,31,30,31,31,30,31,30,31];
% Create two lists of cumulative days for regular and leap years
cum_days = cumsum([0,days]);
cum_days_leap = cumsum([0,days_leap]);
% Check if the start date belongs to a leap year
ind_leap = leapyear(year);
% If the start date belongs to a regular year...
if ind_leap==0
    % The day of the year (doy) is equal to the total number of cumulative
    % days corresponding to the given month for a regular year, plus the
    % day of the starting date
    doy = cum_days(month) + day;
else
    % The day of the year (doy) is equal to the total number of cumulative
    % days corresponding to the given month for a leap year, plus the day
    % of the starting date
    doy = cum_days_leap(month) + day;
end
% Initialize the output date variable
date = zeros(daynum,6);
% For every day that must go in the date list...
for ii=1:daynum
    % Check if the year is regular or leap
    ind_leap = leapyear(year);
    % If the date belongs to a regular year...
    if ind_leap == 0
        % Calculate the month of the corresponding day of the year (doy) as
        % the total number of months (12) minus the number of months ahead
        % of the considered doy, plus 1 because we start counting the
        % months from 1 (January). Uses cum_days for regular years.
        month = 12-sum(doy<=cum_days)+1;
        % Calculate the day subtracting the cumulative number of days of
        % the month calculated at previous step from the day of the year
        day = doy-cum_days(month);
        % Create the date, using the same hours, minutes, and seconds as
        % the starting date
        date(ii,:) = [year,month,day,hours,minutes,seconds];
        % If the day of the year is 365...
        if doy==365
            % Move to the next year
            year = year+1;
            % And reset the day of the year
            doy = 1;
        % ... otherwise
        else
            % Move to the next day of the year
            doy = doy+1;
        end
    else
        % Calculate the month of the corresponding day of the year (doy) as
        % the total number of months (12) minus the number of months ahead
        % of the considered doy, plus 1 because we start counting the
        % months from 1 (January). Uses cum_days_leap for leap years.
        month = 12-sum(doy<=cum_days_leap)+1;
        % Calculate the day subtracting the cumulative number of days of
        % the month calculated at previous step from the day of the year
        day = doy-cum_days_leap(month);
        % Create the date, using the same hours, minutes, and seconds as
        % the starting date
        date(ii,:) = [year,month,day,hours,minutes,seconds];
        % If the day of the year is 366...
        if doy==366
            % Move to the next year
            year = year+1;
            % And reset the day of the year
            doy = 1;
        % ... otherwise
        else
            % Move to the next day of the year
            doy = doy+1;
        end
    end
end
