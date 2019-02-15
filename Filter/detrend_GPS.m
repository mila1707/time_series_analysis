function GPS_dC = detrend_GPS(GPS_DecYr,GPS_C)
% detrend the time series
% using the polyfit to fit the time series
x=GPS_DecYr;
P = polyfit(x,GPS_C,1);
GPS_dC = GPS_C - polyval(P,x);
end

