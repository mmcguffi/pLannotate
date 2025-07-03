"""
Detect fragment features in pLannotate hits
"""
import pandas as pd


def is_fragment(row):
    """Determine if a feature is a fragment"""
    if row["Type"] == "CDS":
        # Full CDS if 100% match
        if row["pi_permatch"] == 100:
            return False
        # Full CDS if length is multiple of 3 and high match percentage
        elif ((row["length"] % 3) == 0) & (row["percmatch"] > 95):
            return False
        else:
            return True
    elif row["Type"] != "CDS":
        # Non-CDS features are fragments if low match percentage
        if row["percmatch"] < 95:
            return True
        else:
            return False
    else:
        # Should not reach here
        raise ValueError(f"Unknown feature type for fragment detection: {row['Type']}")


def detect_fragments(df):
    """Add fragment column to dataframe"""
    if df.empty:
        return df
    
    # Apply fragment detection
    df["fragment"] = df.apply(is_fragment, axis=1)
    
    return df


def main(input_file, output_file):
    """Main function to detect fragments"""
    # Read cleaned hits
    df = pd.read_csv(input_file, sep="\t")
    
    # Detect fragments
    df_with_fragments = detect_fragments(df)
    
    # Save results
    df_with_fragments.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    main(
        input_file=snakemake.input.cleaned,
        output_file=snakemake.output[0]
    )