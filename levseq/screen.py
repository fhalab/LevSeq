import pandas as pd
import numpy as np

# Load the data
file_path = 'path_to_your_csv_file.csv'
data = pd.read_csv(file_path)

# Filter out rows with '*' or '-' in the "Substitutions" column
filtered_data = data[~data["Substitutions"].str.contains(r"[*-]", na=False)]

# Rename columns
filtered_data = filtered_data.rename(columns={"Substitutions": "Sample name", "Well": "Vial"})

# Group by plate and randomly select 8 'PARENT' rows per plate
parent_data = filtered_data[filtered_data['Sample name'] == '#PARENT#']
filtered_parent_data = parent_data.groupby('Plate').apply(lambda x: x.sample(min(8, len(x)))).reset_index(drop=True)

# Update the "Sample name" column to include plate information
filtered_parent_data['Sample name'] = filtered_parent_data['Plate'] + '_' + filtered_parent_data['Sample name']

# Combine the filtered parent data with the rest of the data
non_parent_data = filtered_data[filtered_data['Sample name'] != '#PARENT#']
non_parent_data['Sample name'] = non_parent_data['Plate'] + '_' + non_parent_data['Sample name']

final_data = pd.concat([filtered_parent_data, non_parent_data])

# Add the new columns
final_data['Action'] = 'Inject'
final_data['Sample type'] = 'Sample'
final_data['Injection source'] = 'HipAls'

# Select relevant columns
final_result = final_data[['Sample name', 'Vial', 'Action', 'Sample type', 'Injection source']]

# Save the resulting DataFrame to a new CSV file
output_path = 'path_to_output_csv_file.csv'
final_result.to_csv(output_path, index=False)

