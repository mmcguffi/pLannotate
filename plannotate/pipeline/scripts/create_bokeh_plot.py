"""
Create interactive Bokeh plot for pLannotate annotations
"""
import pandas as pd
from pathlib import Path
import sys

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent.parent))
from bokeh_plot import get_bokeh
from bokeh.embed import file_html
from bokeh.resources import CDN


def main(annotations_file, output_file, is_linear):
    """Main function to create Bokeh plot"""
    # Read annotations
    df = pd.read_csv(annotations_file, sep="\t")
    
    # Create Bokeh plot
    plot = get_bokeh(df, linear=is_linear)
    
    # Generate HTML
    html = file_html(plot, CDN, "pLannotate Plasmid Map")
    
    # Save HTML file
    with open(output_file, "w") as f:
        f.write(html)


if __name__ == "__main__":
    main(
        annotations_file=snakemake.input.annotations,
        output_file=snakemake.output[0],
        is_linear=snakemake.params.linear
    )