%% Process weather data
% Downloaded from https://nsrdb.nrel.gov/data-viewer, Meteostat Prime
% Meridian Europe + Africa dataset, GHI and windspeed, 2005 - 2022, 
% 60 minute intervals 

years = [2005:2022];
windspeeds = cell(12, 24);
irradiance = cell(12, 24);

for year = years
    data = readtable("51.93_4.5/rotterdam_" + string(year) + ".csv");
    for i = 1:height(data)
        month = table2array(data(i,2));
        hour = table2array(data(i,4)) + 1;
        windspeeds{month, hour} = [windspeeds{month, hour}; table2array(data(i, 7));];
        irradiance{month, hour} = [irradiance{month, hour}; table2array(data(i, 6))/1000;];
    end
end

%% Fit the data to probability distributions and save to .mat file

weibull_params_a = cell(12, 24);
weibull_params_b = cell(12, 24);
beta_params_a = cell(12, 24);
beta_params_b = cell(12, 24);
beta_params_scaling = cell(12,24);

for month = 1:12
    for hour = 1:24
        data_wind = windspeeds{month, hour};
        weibull_params = fitdist(data_wind, 'Weibull');
        weibull_params_a{month, hour} = weibull_params.a;
        weibull_params_b{month, hour} = weibull_params.b;

        data_irradiance = irradiance{month, hour};
        beta_params_scaling{month, hour} = max(data_irradiance);
        % 0 solar irradiance in the middle of the night
        if (beta_params_scaling{month, hour} == 0)
            beta_params_a{month, hour} = 0;
            beta_params_b{month, hour} = 0;
            continue
        end
        beta_params = fitdist(data_irradiance./max(data_irradiance), 'Beta');
        beta_params_a{month, hour} = beta_params.a;
        beta_params_b{month, hour} = beta_params.b;
    end
end

%% Save to files for later use
weibull_filename = '51.93_4.5/weibull_params.mat';
save(weibull_filename, 'weibull_params_a', 'weibull_params_b');
beta_filename = '51.93_4.5/beta_params.mat';
save(beta_filename, 'beta_params_a', 'beta_params_b', "beta_params_scaling");


%% Plotting wind speed and irradiance levels for given month/hour, along with fitted distribution

month = 7;
hour = 12;
bins = 20;

data_wind = windspeeds{month, hour};
x_values = [0:0.1:max(data_wind)+5];
y_values = wblpdf(x_values, weibull_params_a{month, hour}, weibull_params_b{month, hour});
bin_width = (max(data_wind) - min(data_wind)) / bins;  % Calculate the bin width
y_values = y_values.*length(data_wind)*bin_width;

figure (1);
histogram(windspeeds{month, hour}, bins);
hold on 
plot(x_values, y_values)
title('Histogram of Wind Speed');
xlabel('Wind Speed (m/s)');
ylabel('Frequency');

data_irradiance = irradiance{month, hour};
x_values = [0.001:0.001:.999];
y_values = betapdf(x_values, beta_params_a{month, hour}, beta_params_b{month, hour});
x_values = x_values.*beta_params_scaling{month, hour};
bin_width = (max(data_irradiance) - min(data_irradiance)) / bins; 
y_values = y_values*bin_width;

figure (2);
histogram(data_irradiance, bins);
hold on 
plot(x_values, y_values)
hold off
title('Histogram of Irradiance');
xlabel('Irradiance (kW/m^2)');
ylabel('Frequency');