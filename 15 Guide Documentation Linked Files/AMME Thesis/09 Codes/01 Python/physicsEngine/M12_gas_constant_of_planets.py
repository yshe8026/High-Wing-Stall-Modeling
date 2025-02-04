# Gas constant of planets in the solar system + Titan
# https://pds-atmospheres.nmsu.edu/education_and_outreach/encyclopedia/gas_constant.htm


# Gas Constant

# Runiv = cp - cv = 8.3143 Joules/mole/K = 8.3143 · 107 erg/mole/K
# <mw>air = 28.964 gm/mole
# Rg(air) = Runiv/<mw>air = 287.05 Joule/Kg/K
# Cp(air) = 1005 J/Kg/K
# kappa ident Rg/Cp = Runiv/cp
# gamma ident Cp/Cv = cp/cv

#           <mw>    Rgas     Cp     kappa   gamma
#           gm/mole J/gm/K  J/gm/K

# Venus	    44.01	0.18892	0.8501	0.2222	1.2857
# Earth	    28.96	0.28710	1.0040	0.2860	1.4005
# Mars	    44.01	0.18892	0.8312	0.2273	1.2941
# Jupiter	2.22	3.74518	12.3591	0.3030	1.4348
# Saturn	2.14	3.89246	14.0129	0.2778	1.3846
# Uranus	2.30	3.61491	13.0137	0.2778	1.3846
# Neptune	2.30	3.61491	13.0137	0.2778	1.3846
# Titan	    28.67	0.29000	1.0440	0.2778	1.3846

import pandas as pd

# Define the data as a dictionary
data = {
    "Planet": ["Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Titan"],
    "mw (gm/mole)": [44.01, 28.96, 44.01, 2.22, 2.14, 2.30, 2.30, 28.67],
    "Rgas (J/gm/K)": [0.18892, 0.28710, 0.18892, 3.74518, 3.89246, 3.61491, 3.61491, 0.29000],
    "Cp (J/gm/K)": [0.8501, 1.0040, 0.8312, 12.3591, 14.0129, 13.0137, 13.0137, 1.0440],
    "kappa": [0.2222, 0.2860, 0.2273, 0.3030, 0.2778, 0.2778, 0.2778, 0.2778],
    "gamma": [1.2857, 1.4005, 1.2941, 1.4348, 1.3846, 1.3846, 1.3846, 1.3846]
}

# Create a DataFrame
df = pd.DataFrame(data)

# Display the DataFrame
print(df)

# Look up the values for Mars
mars_data = df[df["Planet"] == "Mars"]

# Extract the Rgas and gamma values
mars_Rgas = mars_data["Rgas (J/gm/K)"].values[0]
mars_gamma = mars_data["gamma"].values[0]

print(f"Rgas for Mars: {mars_Rgas} J/gm/K")
print(f"Gamma for Mars: {mars_gamma}")