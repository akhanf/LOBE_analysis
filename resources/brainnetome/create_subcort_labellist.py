import csv

# Input and output file paths
csv_file = "subregion_func_network_Yeo_updated.csv"
output_file = "subcort_labellist.txt"

# Process the input CSV file
with open(csv_file, 'r') as infile, open(output_file, 'w') as outfile:
    reader = csv.reader(infile)
    for row in reader:
        # Extract relevant data from columns
        index = row[0]
        name = row[1]
        detailed_name = row[2]
        
        # Determine L_ or R_ based on the presence of "L" or "R" in the third column
        if "L" in detailed_name:
            side = "L_"
        elif "R" in detailed_name:
            side = "R_"
        else:
            raise ValueError(f"Invalid side information in detailed_name: {detailed_name}")
        
        # Construct the new name
        new_name = f"{side}{name}_Subcort"
        
        # Write the output in the required format
        outfile.write(new_name + "\n")
        outfile.write(f"{index} 128 128 128 255\n")  # Example values, can adjust as needed

print(f"Updated file saved to {output_file}")


