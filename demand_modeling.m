% Downloaded from https://data.open-power-system-data.org/household_data/

opts = detectImportOptions("household_data_60min_singleindex_filtered-2.csv");
data = readtable("household_data_60min_singleindex_filtered-2.csv", 'ReadVariableNames', true, 'TextType', 'string', 'Delimiter', ',');
header_names = data.Properties.VariableNames

processed_data = table();

for i = 4:length(header_names)
    header_name = string(header_names(i))
    diff_value = [NaN; diff(data.(header_name))];
    processed_data.(header_name) = diff_value;
end
    
date = datetime(data.("utc_timestamp"), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ssZ', 'TimeZone','UTC');
processed_data.("year") = year(date);
processed_data.("month") = month(date);
processed_data.("day") = day(date);
processed_data.("hour") = hour(date);


years = get_valid_years(processed_data, "DE_KN_residential4_grid_import");
date = datetime("2016-06-28", "Format", "yyyy-MM-dd");
isweekend(date)
demand_industrial = get_energy_demand_day(processed_data, 2016, 6, 23, "DE_KN_industrial1_grid_import");
demand_public1 = get_energy_demand_day(processed_data, 2016, 6, 23, "DE_KN_public1_grid_import");
demand_residential1 = get_energy_demand_day(processed_data, 2016, 6, 23, "DE_KN_residential1_grid_import");
demand_residential2 = get_energy_demand_day(processed_data, 2016, 6, 23, "DE_KN_residential2_grid_import");
demand_residential3 = get_energy_demand_day(processed_data, 2016, 6, 23, "DE_KN_residential3_grid_import");
demand_residential4 = get_energy_demand_day(processed_data, 2016, 6, 23, "DE_KN_residential4_grid_import");
demand_residential5 = get_energy_demand_day(processed_data, 2016, 6, 23, "DE_KN_residential5_grid_import");
demand_residential6 = get_energy_demand_day(processed_data, 2016, 6, 23, "DE_KN_residential6_grid_import");

function years = get_valid_years(processed_data, param_name)
    uniqueYears = unique(processed_data.("year"));
    years = [];
    for i = 1:length(uniqueYears)
        curr_year = uniqueYears(i);    
        year_data = processed_data(processed_data.("year") == curr_year, :);
        if all(~isnan(year_data.(param_name))) 
            years = [years; curr_year]; 
        end
    end
end

function demand = get_energy_demand_day(processed_data, year, months, days, param_name)
    demand = [];
    year_demand = processed_data(processed_data.("year") == year, :);
    for month = months
        month_demand = year_demand(year_demand.("month") == month, :);
        for day = days
            day_demand = month_demand(month_demand.("day") == day, :);
            demand = [demand; day_demand.(param_name)];
        end
    end
end

figure(1)
plot(demand_industrial)
hold on
plot(demand_public1)
legend(["Industrial1", "Public1"])
title("Industrial power demand")

figure(2)
% plot(demand_residential1)
% hold on
% plot(demand_residential2)
% hold on
% plot(demand_residential4)
% hold on
% plot(demand_residential5)
% hold on 
plot(demand_residential1+demand_residential2+demand_residential3+demand_residential4+demand_residential5+demand_residential6)
xlabel("Hour")
ylabel("Power demand (kW)")
title("Residential power demand")
