#!/usr/bin/env python3
"""
Create databases.yml configuration file from version information using a template.

Usage:
    python3 create_databases_yml.py \
        --rfam-version rfam_version.txt \
        --snapgene-version snapgene_version.txt \
        --fpbase-version fpbase_version.txt \
        --swissprot-version swissprot_version.txt \
        --output databases.yml \
        --template templates/databases_template.yml
"""

import argparse
from datetime import datetime
from pathlib import Path


def read_version_file(version_file: Path) -> str:
    """Read version information from a file."""
    try:
        with open(version_file, "r") as f:
            return f.read().strip()
    except Exception:
        return f"Unknown version ({datetime.now().strftime('%Y-%m-%d')})"


def create_databases_yml_from_template(
    template_file: Path,
    rfam_version: str,
    snapgene_version: str,
    fpbase_version: str,
    swissprot_version: str,
) -> str:
    """Create the databases.yml content from template file."""
    try:
        with open(template_file, "r") as f:
            template_content = f.read()

        # Replace version placeholders
        yml_content = template_content.format(
            rfam_version=rfam_version,
            snapgene_version=snapgene_version,
            fpbase_version=fpbase_version,
            swissprot_version=swissprot_version,
        )

        return yml_content

    except FileNotFoundError:
        raise FileNotFoundError(f"Template file not found: {template_file}")
    except Exception as e:
        raise Exception(f"Error processing template: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Create databases.yml configuration file"
    )
    parser.add_argument("--rfam-version", required=True, help="Rfam version file")
    parser.add_argument(
        "--snapgene-version", required=True, help="SnapGene version file"
    )
    parser.add_argument("--fpbase-version", required=True, help="FPbase version file")
    parser.add_argument(
        "--swissprot-version", required=True, help="Swiss-Prot version file"
    )
    parser.add_argument("--output", required=True, help="Output databases.yml file")
    parser.add_argument(
        "--template",
        default="templates/databases_template.yml",
        help="Template YAML file (default: templates/databases_template.yml)",
    )

    args = parser.parse_args()

    # Read version information
    rfam_version = read_version_file(Path(args.rfam_version))
    snapgene_version = read_version_file(Path(args.snapgene_version))
    fpbase_version = read_version_file(Path(args.fpbase_version))
    swissprot_version = read_version_file(Path(args.swissprot_version))

    # Create YAML content from template
    template_file = Path(args.template)
    yml_content = create_databases_yml_from_template(
        template_file, rfam_version, snapgene_version, fpbase_version, swissprot_version
    )

    # Write output file
    output_file = Path(args.output)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w") as f:
        f.write(yml_content)

    print(f"Created databases.yml: {output_file}")
    print(f"  Rfam: {rfam_version}")
    print(f"  SnapGene: {snapgene_version}")
    print(f"  FPbase: {fpbase_version}")
    print(f"  Swiss-Prot: {swissprot_version}")


if __name__ == "__main__":
    main()
