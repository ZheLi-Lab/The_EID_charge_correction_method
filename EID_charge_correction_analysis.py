import os
import pandas as pd
import numpy as np

def get_sorted_csv_files(directory, prefix='state_s', suffix='.csv'):
    """Retrieve CSV files sorted by numeric suffix."""
    csv_files = [f for f in os.listdir(directory) if f.startswith(prefix) and f.endswith(suffix)]
    csv_files.sort(key=lambda x: int(x.split('_')[1].split('.')[0][1:]))
    return csv_files

def process_csv_file(file_path):
    """Process a single CSV file, extracting data and lambda values."""
    df = pd.read_csv(file_path, sep='|')
    columns = df.columns[4:]
    lambda_list = [col.split(',')[1].replace(')', '').strip() for col in columns]
    data = df.iloc[:, 4:].values  # Extract data starting from the 6th column
    return df, data, list(map(float, lambda_list))

def solve_for_ABC(lambdas, Us):
    """Fit a quadratic equation U = A*lambda^2 + B*lambda + C and compute RMSE."""
    lambdas = np.array(lambdas)
    Us = np.array(Us)
    X = np.vstack((lambdas**2, lambdas, np.ones_like(lambdas))).T
    coeffs, _, _, _ = np.linalg.lstsq(X, Us, rcond=None)
    A, B, C = coeffs
    Us_fit = X @ coeffs  # Compute fitted values
    rmse = np.sqrt(np.mean((Us - Us_fit)**2))  # Calculate RMSE
    return A, B, C, rmse

def process_all_files(directory):
    """Process all files in the directory, fitting data and saving results."""
    csv_files = get_sorted_csv_files(directory)

    for csv_file in csv_files:
        file_path = os.path.join(directory, csv_file)
        print(f"Processing file: {csv_file}")

        # Process the CSV file
        original_df, data, lambda_list = process_csv_file(file_path)

        # Calculate A, B, C, and RMSE for each row
        abc_results = []
        for row in data:
            A, B, C, rmse = solve_for_ABC(lambda_list, row)
            equation = f"U = {A:.4f}*lambda^2 + {B:.4f}*lambda + {C:.4f}"
            abc_results.append({'A': A, 'B': B, 'C': C, 'Equation': equation, 'RMSE': rmse})

        # Save ABC results to a CSV file
        abc_df = pd.DataFrame(abc_results)
        abc_file = csv_file.replace('state_s', 'ABC_state_s')
        abc_df.to_csv(abc_file, index=False)
        print(f"ABC results saved to {abc_file}")

        # Calculate Delta U
        base_lambda = lambda_list[0]
        delta_lambda = [l - base_lambda for l in lambda_list[1:]]
        #delta_lambda = [lambda_list[i] - lambda_list[i - 1] for i in range(1, len(lambda_list))]
        B_values = abc_df['B'].values
        delta_us = np.zeros((len(data), len(delta_lambda)))

        for i, b in enumerate(B_values):
            delta_us[i, :] = b * np.array(delta_lambda)

        # Create Delta U DataFrame
        delta_us_file = csv_file.replace('state_s', 'state_ΔUs')
        delta_us_df = original_df.copy()  # Start with the original DataFrame

        # Ensure at least 6 columns exist, add a zero-filled column for the 6th column
        if len(delta_us_df.columns) < 5:
            delta_us_df.insert(4, "Delta_Zero", 0)
        else:
            delta_us_df.iloc[:, 4] = 0

        # Fill calculated Delta U values starting from the 7th column
        for i in range(len(delta_lambda)):
            delta_us_df[delta_us_df.columns[5 + i]] = delta_us[:, i]

        delta_us_df.to_csv(delta_us_file, sep='|', index=False)
        print(f"ΔU results saved to {delta_us_file}")

if __name__ == "__main__":
    directory = './'  # Update with the directory containing your CSV files
    process_all_files(directory)

