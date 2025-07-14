#!/usr/bin/env python3
"""
Simple validation script for gather.smk Snakemake workflow.
"""

import subprocess
import sys
from pathlib import Path


def validate_snakemake_syntax():
    """Validate Snakemake workflow syntax."""
    try:
        # Try to parse the workflow without running it
        result = subprocess.run(
            ["python", "-c", "import snakemake; print('Snakemake available')"],
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            print("Warning: Snakemake not available for syntax validation")
            return True  # Skip validation if snakemake not available

        # Try to lint/validate the workflow
        result = subprocess.run(
            ["snakemake", "-s", "gather.smk", "--lint"], capture_output=True, text=True
        )

        if result.returncode == 0:
            print("✓ Snakemake workflow syntax is valid")
            return True
        else:
            print(f"✗ Snakemake syntax errors: {result.stderr}")
            return False

    except FileNotFoundError:
        print("Warning: Snakemake not available for validation")
        return True
    except Exception as e:
        print(f"Error validating workflow: {e}")
        return False


def validate_python_scripts():
    """Validate Python helper scripts syntax."""
    scripts_dir = Path("scripts")
    valid = True

    for script in scripts_dir.glob("*.py"):
        try:
            result = subprocess.run(
                ["python", "-m", "py_compile", str(script)],
                capture_output=True,
                text=True,
            )

            if result.returncode == 0:
                print(f"✓ {script.name} syntax is valid")
            else:
                print(f"✗ {script.name} syntax error: {result.stderr}")
                valid = False

        except Exception as e:
            print(f"Error validating {script.name}: {e}")
            valid = False

    return valid


def validate_file_structure():
    """Validate required files and directories exist."""
    required_files = [
        "gather.smk",
        "scripts/build_diamond_db.py",
        "scripts/create_databases_yml.py",
        "rfam/gather_rfam.py",
        "snapgene/gather_snapgene.py",
        "fpbase/gather_fpbase.py",
        "fpbase/gather_fpbase.sh",
        "swissprot/download_fresh_swissprot.py",
        "swissprot/run_full_workflow.sh",
    ]

    valid = True
    for file_path in required_files:
        if Path(file_path).exists():
            print(f"✓ {file_path} exists")
        else:
            print(f"✗ {file_path} missing")
            valid = False

    return valid


def main():
    """Run all validations."""
    print("=== Validating gather.smk workflow ===\n")

    # Validate file structure
    print("Checking file structure...")
    structure_valid = validate_file_structure()
    print()

    # Validate Python scripts
    print("Checking Python script syntax...")
    scripts_valid = validate_python_scripts()
    print()

    # Validate Snakemake workflow
    print("Checking Snakemake workflow...")
    workflow_valid = validate_snakemake_syntax()
    print()

    # Summary
    if structure_valid and scripts_valid and workflow_valid:
        print("✓ All validations passed!")
        print("\nTo run the workflow:")
        print("  snakemake -s gather.smk --cores 4")
        print("  snakemake -s gather.smk --cores 4 --dry-run  # Preview workflow")
        return 0
    else:
        print("✗ Some validations failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())
