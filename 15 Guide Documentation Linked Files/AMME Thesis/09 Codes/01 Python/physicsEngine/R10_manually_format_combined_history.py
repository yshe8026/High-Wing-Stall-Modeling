import pandas as pd

# Specify the correct file path
file_path = 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\richard_checked_terms\\combined_history.txt'

# Load the data
df = pd.read_csv(file_path, sep=r'\s+')

# Drop duplicate rows based on 'global_time[s]', keeping the last occurrence
df_cleaned = df.drop_duplicates(subset=['global_time[s]'], keep='last')

# Specify the output file path
output_path = 'C:\\Users\\Richard\\00 Richard Apps\\Steam\\steamapps\\common\\X-Plane 11\\richard_checked_terms\\formatted_combined_history.txt'

# Define the column width
cwidth = 25

# Function to format the columns
def format_line(global_time, global_dt, dX_flow_sep_by_dt, term_0, X_flow_sep, term_1, term_2, term_3, term_4, term_5, alpha, alpha_dot, a_1, alpha_star, tau_1, tau_2, Vt, cwidth):
    # Format global_time to 5 significant figures
    formatted_global_time = f"{global_time:.5g}"
    return f"{formatted_global_time:<{cwidth}}{global_dt:<{cwidth}}{dX_flow_sep_by_dt:<{cwidth}}{term_0:<{cwidth}}{X_flow_sep:<{cwidth}}{term_1:<{cwidth}}{term_2:<{cwidth}}{term_3:<{cwidth}}{term_4:<{cwidth}}{term_5:<{cwidth}}{alpha:<{cwidth}}{alpha_dot:<{cwidth}}{a_1:<{cwidth}}{alpha_star:<{cwidth}}{tau_1:<{cwidth}}{tau_2:<{cwidth}}{Vt:<{cwidth}}\n"

# Write the cleaned data to the output file with aligned columns
with open(output_path, 'w') as file:
    # Write the column titles
    file.write(f"{'global_time[s]':<{cwidth}}{'global_dt[s]':<{cwidth}}{'dX_flow_sep_by_dt[-/s]':<{cwidth}}{'term_0[-/s]':<{cwidth}}{'X_flow_sep[-]':<{cwidth}}{'term_1[-/s]':<{cwidth}}{'term_2[-/s]':<{cwidth}}{'term_3[-]':<{cwidth}}{'term_4[-]':<{cwidth}}{'term_5[-]':<{cwidth}}{'alpha[rad]':<{cwidth}}{'alpha_dot[rad/s]':<{cwidth}}{'a_1[-]':<{cwidth}}{'alpha_star[rad]':<{cwidth}}{'tau_1[s]':<{cwidth}}{'tau_2[s]':<{cwidth}}{'Vt[m/s]':<{cwidth}}\n")
    for index, row in df_cleaned.iterrows():
        file.write(format_line(row['global_time[s]'], row['global_dt[s]'], row['dX_flow_sep_by_dt[-/s]'], row['term_0[-/s]'], row['X_flow_sep[-]'], row['term_1[-/s]'], row['term_2[-/s]'], row['term_3[-]'], row['term_4[-]'], row['term_5[-]'], row['alpha[rad]'], row['alpha_dot[rad/s]'], row['a_1[-]'], row['alpha_star[rad]'], row['tau_1[s]'], row['tau_2[s]'], row['Vt[m/s]'], cwidth))

print(f"Formatted data saved to {output_path}")