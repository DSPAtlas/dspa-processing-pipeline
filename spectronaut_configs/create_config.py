import pandas as pd

def generate_config(names, replicate_number):
    base_output_path = "G:\\Biognosys\\Spectronaut\\results\\"
    settings_file = "G:/Biognosys/1156159636_BatchResources/dynaprot_directdia_settings.prop"
    condition_file =  "G:/Biognosys/1156159636_BatchResources/conditions.tsv"
    fasta_file = "g:/biognosys/spectronaut/fasta/human9606241108.bgsfasta"
    raw_files_dir =  "Y:\\LiP-Atlas\\data\\PXD015446_Ilaria\\Raw_data\\"
    
    # Open the configurations file to write configuration details
    with open("configurations.txt", 'w') as config_file:
        # Also open the conditions file to write conditions details
        with open("conditions.csv", 'w') as conditions_file:
            # Write header for conditions.csv
            conditions_file.write("Reference;Run Label;Condition;Fraction;Replicate;Quantity Correction Factor;Label;Color;File Name\n")
            
            # Process each name
            for idx, name in enumerate(names):
               
                part_number = replicate_number[idx]  

                # Write to configuration file

                config_file.write(f' -direct\n')
                config_file.write(f'-s "{settings_file}"\n')
                config_file.write(f'-con "{condition_file}"\n')
                config_file.write(f'-n "PXD015446_Ilaria_FK506_{part_number}"\n')
                config_file.write(f'-o "{base_output_path}PXD015446_Ilaria_FK506"\n')
                config_file.write(f'-fasta "{fasta_file}"\n')
                config_file.write(f'-r "{raw_files_dir}{name}.raw"\n\n')
                
                # Prepare and write a line to the conditions file
                conditions_line = f"False;{name};{part_number};NA;1;1;Not Defined;#505050;{name}\n"
                conditions_file.write(conditions_line)

# Example usage of the function
df = pd.read_csv("experiments.csv", sep=";", encoding='iso-8859-1')
experiment_name = "PXD015446_Ilaria_FK506"
filtered_df = df[df['Experiment Name'] == experiment_name]
raw_file_names = filtered_df['Raw File Name'].tolist()
replicate_number = filtered_df['Replicate Number'].tolist()
generate_config(raw_file_names, replicate_number)
