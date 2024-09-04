
 clear Place;
 % The amplitude H and phase for each constituent are determined from the tidal record by least-squares
 % fitting to the observations of the amplitudes of the astronomical terms with expected frequencies and phases.
 % The number of constituents needed for accurate prediction varies from place to place.
 % In making up the tide tables for Long Island Sound, the National Oceanic and Atmospheric Administration
 % uses 23 constituents. The eleven whose amplitude is greater than .1 foot are:
 Place(1).Name = 'Bridgeport, Cn';	% Counting time in hours from midnight starting Sunday 1 September 1991.
 %                  M2      S2     N2    L2    K2     K1     O1    Sa   nu2    Mm     P1...
 Place(1).A = [  3.185   0.538  0.696 0.277 0.144  0.295  0.212 0.192 0.159 0.108  0.102];	% Tidal heights (feet)
 Place(1).P = [-127.24 -343.66 263.60 -4.72 -2.55 142.02 505.93 301.5 45.70 86.82 340.11];	% Phase (degrees). 
 % The values for these coefficients are taken from http://www.math.sunysb.edu/~tony/tides/harmonic.html
 % which originally came from a table published by the US. National Oceanic and Atmospheric Administration.

 % Calculate a tidal height curve, in terms of hours since the start time.
 PlaceCount = 1;
 Colour=cellstr(strvcat('g','r','b','c','m','y','k'));	% A collection.
 clear y;
 step = 0.125; LastHour = 720; % 8760 hours in a year.
 n = LastHour/step + 1;
 y(1:n,1:PlaceCount) = 0;
 t = (0:step:LastHour)/24;
 for it = 1:PlaceCount
  i = 0;
  for h = 0:step:LastHour
   i = i + 1;
   y(i,it) = sum(Place(it).A.*cosd(Tide.Speed*h + Place(it).P)); %Sum terms A(j)*cos(speed(j)*h + p(j)) for j = 1:Tide.Constituents.
  end;      % Should use cos(ix) = 2*cos([i - 1]*x)*cos(x) - cos([i - 2]*x), but, for clarity...
 end;

 figure(1); clf; hold on; title('Tidal Height'); xlabel('Days');
 for it = 1:PlaceCount
  plot(t,y(1:n,it),Colour{it});
 end;
 legend(Place(1:PlaceCount).Name,'Location','NorthWest');

