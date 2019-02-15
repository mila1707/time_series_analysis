function [GPS_Ni2,GPS_Ei2,GPS_Ui2,decdatei2] = GPS_interpolation(GPS_DecYr,GPS_N,GPS_E,GPS_U,range1,range2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
year2=[2008*ones(yeardays(2008),1);
       2009*ones(yeardays(2009),1);
       2010*ones(yeardays(2010),1);
       2011*ones(yeardays(2011),1);
       2012*ones(yeardays(2012),1);
       2013*ones(yeardays(2013),1);
       2014*ones(yeardays(2014),1);
       2015*ones(yeardays(2015),1);
       2016*ones(yeardays(2016),1);
       2017*ones(yeardays(2017),1);];
%        2018*ones(yeardays(2018),1);];
doy2=[[1:yeardays(2008)]';
      [1:yeardays(2009)]';
      [1:yeardays(2010)]';
      [1:yeardays(2011)]';
      [1:yeardays(2012)]';
      [1:yeardays(2013)]';
      [1:yeardays(2014)]';
      [1:yeardays(2015)]';
      [1:yeardays(2016)]';
      [1:yeardays(2017)]';];
%       [1:yeardays(2018)]';];
decdatei2=year2+doy2./yeardays(year2);

GPS_Ni2 = interp1(GPS_DecYr,GPS_N,decdatei2,'linear');
GPS_Ei2 = interp1(GPS_DecYr,GPS_E,decdatei2,'linear');
GPS_Ui2 = interp1(GPS_DecYr,GPS_U,decdatei2,'linear');

end

