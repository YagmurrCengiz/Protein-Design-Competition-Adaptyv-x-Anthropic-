import pandas as pd

input_path = (
    "data/structures/"
    "egfr_high_priority_analysis.csv"
)

output_path = (
    "data/structures/"
    "egfr_candidate_domains.csv"
)

df = pd.read_csv(input_path)

# These are the commonly used approximate boundaries
# for the four EGFR extracellular domains.

def assign_domain(position):

    if 25 <= position <= 190:
        return "Domain I"

    elif 191 <= position <= 292:
        return "Domain II"

    elif 293 <= position <= 444:
        return "Domain III"

    elif 445 <= position <= 645:
        return "Domain IV"

    else:
        return "Outside extracellular region"

df["EGFR_domain"] = df["position"].apply(
    assign_domain
)

print()
print("=" * 90)
print("EGFR CANDIDATE DOMAIN ANALYSIS")
print("=" * 90)

print()

for _, row in df.iterrows():

    print(
        f"{int(row['position'])}: "
        f"{row['human']}->{row['mouse']} | "
        f"{row['EGFR_domain']}"
    )

df.to_csv(
    output_path,
    index = False
)

print()
print("=" * 90)
print("DOMAIN ANALYSIS SAVED")
print("=" * 90)

print()

print(
    "Output:",
    output_path
)

print(
    "Candidates:",
    len(df)
)