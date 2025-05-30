function cvc = contantsvalues_convertions()
    
    global radians_to_degrees radians_to_degrees Earth_gravitational_constant Earth_radius

    cvc.degrees_to_radians           = pi / 180.0;
    cvc.radians_to_degrees           = 180.0 / pi;
    cvc.Earth_gravitational_constant = (6.67430*10^-11) * ((5.9722+0.0006)*10^24); 
    cvc.Earth_radius                 = ((2 * 6378.1370) + 6356.7523)/3;
    
end
