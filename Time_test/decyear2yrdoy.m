function [yr,doy] = decyear2yrdoy(timeline)

% decyear2yrdoy.m converts a timeline given in decimal year into two
% variables containing the year and the corresponding day of the year
% (doy).
% -------------------------------------------------------------------------
% INPUT
% timeline: Vector containing the date in decimal year. Size: Tx1
% -------------------------------------------------------------------------
% OUTPUT
% yr: Vector containing the year. Size: Tx1
% doy: Vector containing the day of the year. Size: Tx1
% -------------------------------------------------------------------------
% 
% Adriano Gualandi - 21 Mar 2017
% California Institute of Technology
% Geological and Planetary Science Division

% Find the year variable
yr = fix(timeline);
% Find the indexes for which the year is a leap year
%ind_leap = leapyear(yr);
ind_leap = leapyear_notlicensed(yr);
% Calculate the day of the year as the remainder between the decimal time
% and the year, times the number of days in one year
doy = (timeline-yr)*365;
% Correct using 366 days in one year for leap years
doy(ind_leap) = (timeline(ind_leap)-yr(ind_leap))*366;