function date_time = decyear2date(timeline)

% decyear2date.m converts times from decimal year format to date format. It
% takes advantage of the two functions decyear2yrdoy.m and yrdoy2date.m.
% -------------------------------------------------------------------------
% INPUT
% timeline: Vector containing the date in decimal year. Size: Tx1
% -------------------------------------------------------------------------
% OUTPUT
% date_time: Matrix containing the year, month, day, hours, minutes, and
%            seconds. Size: Tx6
% -------------------------------------------------------------------------
% 
% Adriano Gualandi - 21 Mar 2017
% California Institute of Technology
% Geological and Planetary Science Division

% Convert the timeline in year and day of the year (doy)
[year,doy] = decyear2yrdoy(timeline);
% Convert the year and doy in date format
date_time = yrdoy2date(year,doy);