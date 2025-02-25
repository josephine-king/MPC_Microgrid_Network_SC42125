%% Process weather data

years = [2005:2022];
windspeeds = cell(12, 1);
irradiance = cell(12, 24);

for year = years
    data = readtable("rotterdam_weather_data/rotterdam_" + string(year) + ".csv");
    for i = 1:height(data)
        month = table2array(data(i,2));
        hour = table2array(data(i,4)) + 1;
        windspeeds{month} = [windspeeds{month}; table2array(data(i, 7));];
        irradiance{month, hour} = [irradiance{month, hour}; table2array(data(i, 6));];
    end
end

%% Fit the data to probability distributions and save to .mat file

weibull_params_a = cell(12, 1);
weibull_params_b = cell(12, 1);
beta_params_a = cell(12, 24);
beta_params_b = cell(12, 24);
beta_params_scaling = cell(12,24);

for month = 1:12
    data_wind = windspeeds{month};
    weibull_params = fitdist(data_wind, 'Weibull');
    weibull_params_a{month} = weibull_params.a;
    weibull_params_b{month} = weibull_params.b;
    for hour = 1:24
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

% Save to files for later use
weibull_filename = 'rotterdam_weather_data/weibull_params.mat';
save(weibull_filename, 'weibull_params_a', 'weibull_params_b');
beta_filename = 'rotterdam_weather_data/beta_params.mat';
save(beta_filename, 'beta_params_a', 'beta_params_b');


%% Plotting wind speed and irradiance levels for given month/hour, along with fitted distribution

month = 10;
hour = 13;
bins = 50;

data_wind = windspeeds{month};
x_values = [0:0.1:max(data_wind)+5];
y_values = wblpdf(x_values, weibull_params_a{month}, weibull_params_b{month});
bin_width = (max(data_wind) - min(data_wind)) / bins;  % Calculate the bin width
y_values = y_values.*length(data_wind)*bin_width;

figure (1);
histogram(windspeeds{month}, bins);
hold on 
plot(x_values, y_values)
title('Histogram of Wind Speed');
xlabel('Wind Speed (m/s)');
ylabel('Frequency');

data_irradiance = irradiance{month, hour};
x_values = [0:0.001:1];
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