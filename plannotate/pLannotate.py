import argparse
import sys

import click
import yaml

import plannotate.resources as rsc
from plannotate.annotate import annotate


@click.group()
@click.version_option(prog_name=__package__)
def main():
    pass


@main.command("run")
@click.option('--yaml_file', default = rsc.get_yaml_path(), help="path to YAML file.", type=click.Path(exists=True))
def main_run(yaml_file): 
    """Runs pLannotate."""
    if rsc.databases_exist():
        run_plannotate(yaml_file)
    else:
        print("Databases not downloaded. Run 'plannotate setupdb' to download databases.")


@main.command("yaml")
def main_yaml():
    """Prints YAML file to stdout for custom database modification."""
    with open (rsc.get_yaml_path(), 'r') as stream:
        print(yaml.dump(yaml.load(stream, Loader = yaml.SafeLoader), default_flow_style=False))
        

@main.command("setupdb")
def main_setupdb():
    """Downloads databases; required for use of pLannotate."""
        
    if rsc.databases_exist():
        print("Databases already downloaded.")
        print()
        
    else:
       rsc.download_databases()

    print("Run 'plannotate run {arguments}' to launch pLannotate.")
    print("To get a list of available arguments for command line use, run 'plannotate --help'.")
    print("Please also consider citing: https://doi.org/10.1093/nar/gkab374 :)")


@main.command("batch")
@click.option("--input","-i", 
                help=f"location of a FASTA or GBK file")
@click.option("--output","-o", default = f"./",  
                help="location of output folder. DEFAULT: current dir")
@click.option("--file_name","-f", default = "",  
                help="name of output file (do not add extension). DEFAULT: input file name")
@click.option("--suffix","-s", default = "_pLann",  
                help="suffix appended to output files. Use '' for no suffix. DEFAULT: '_pLann'")
@click.option("--yaml_file","-y", default=rsc.get_yaml_path(), 
               help="path to YAML file for custom databases. DEFAULT: builtin")
@click.option("--linear","-l", is_flag=True, 
                help="enables linear DNA annotation")
@click.option("--html","-h", is_flag=True, 
                help="creates an html plasmid map in specified path")
@click.option("--csv","-c", is_flag=True, 
                help="creates a cvs file in specified path")
@click.option("--detailed","-d", is_flag=True, 
                help="uses modified algorithm for a more-detailed search with more false positives")
@click.option("--no_gbk","-x", is_flag=True, 
                help="supresses GenBank output file")
def main_batch(**kwargs):
    """
    Annotates engineered DNA sequences, primarily plasmids. Accepts a FASTA or GenBank file and outputs
    a GenBank file with annotations, as well as an optional interactive plasmid map as an HTLM file.
    """
    if not rsc.databases_exist():
        print("Databases not downloaded. Run 'plannotate setupdb' to download databases.")
        sys.exit()

    name, ext = rsc.get_name_ext(kwargs['input'])

    if kwargs['file_name'] == "":
        kwargs['file_name'] = name

    inSeq = rsc.validate_file(kwargs['input'], ext, max_length = float("inf"))

    recordDf = annotate(inSeq, kwargs['yaml_file'], kwargs['linear'], kwargs['detailed'])

    if kwargs['no_gbk'] == False:
        gbk = rsc.get_gbk(recordDf, inSeq, kwargs['linear'])
        with open(f"{kwargs['output']}/{kwargs['file_name']}{kwargs['suffix']}.gbk", "w") as handle:
            handle.write(gbk)

    if kwargs['csv']:
        csv_df = rsc.get_clean_csv_df(recordDf)
        csv_df.to_csv(f"{kwargs['output']}/{kwargs['file_name']}{kwargs['suffix']}.csv", index = None)


if __name__ == '__main__':
    main()
