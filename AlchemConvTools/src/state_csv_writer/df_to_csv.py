import os
import pandas as pd

def large_df_to_csvs(example_df, output_dir):
    """
    Split a large DataFrame containing multiple states into individual DataFrames for each state and output them as CSV files.

    :param example_df: The large DataFrame containing multiple states.
    :param output_dir: The directory where the CSV files will be saved.
    """
    # Get the list of lambda information
    lambda_info_list = list(example_df.index.names)
    lambda_info_list.pop(0)

    # Group the DataFrame by lambda information
    every_single_csv_list = [j.dropna(axis=1, how='all') for i, j in example_df.groupby(lambda_info_list, sort=False)]

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Iterate over each state's DataFrame and output as CSV file
    for i, df in enumerate(every_single_csv_list, start=1):
        # Convert column names to string format
        df.columns = pd.Index([str(col) for col in df.columns])

        # Generate the CSV file name
        csv_filename = f"state_s{i}.csv"
        csv_path = os.path.join(output_dir, csv_filename)

        # Output the DataFrame as a CSV file
        df.to_csv(csv_path, sep='|', index=True, header=True)

#    print(f"Successfully output {len(every_single_csv_list)} state DataFrames as CSV files.")
