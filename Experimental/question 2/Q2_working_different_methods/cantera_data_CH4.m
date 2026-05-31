function [phi,sf] = cantera_data_CH4()

data = [
    0.8     0.274326
    1.0     0.380927
    1.2     0.336674

];

phi = data(:,1);
sf = data(:,2)*1000;

end