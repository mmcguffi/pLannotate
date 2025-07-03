#!/usr/bin/env python3
"""
Test script to verify the Snakemake pipeline produces
the same results as the original implementation.
"""

import sys
import pandas as pd
from tempfile import TemporaryDirectory

# Import both versions
from annotate import annotate as new_annotate
from annotate_original import annotate as original_annotate


def compare_dataframes(df1, df2, name1="New", name2="Original"):
    """Compare two DataFrames and report differences."""
    print(f"\nComparing {name1} vs {name2}:")
    print(f"  {name1} shape: {df1.shape}")
    print(f"  {name2} shape: {df2.shape}")
    
    if df1.shape != df2.shape:
        print("  ❌ Different shapes!")
        return False
    
    # Compare columns
    if not set(df1.columns) == set(df2.columns):
        print("  ❌ Different columns!")
        print(f"    {name1} only: {set(df1.columns) - set(df2.columns)}")
        print(f"    {name2} only: {set(df2.columns) - set(df1.columns)}")
        return False
    
    # Sort by key columns for comparison
    sort_cols = ['qstart', 'qend', 'sseqid', 'db']
    df1_sorted = df1.sort_values(sort_cols).reset_index(drop=True)
    df2_sorted = df2.sort_values(sort_cols).reset_index(drop=True)
    
    # Compare values
    differences = []
    for col in df1.columns:
        if df1_sorted[col].dtype in ['float64', 'int64']:
            # Numeric comparison with tolerance
            if not pd.np.allclose(df1_sorted[col], df2_sorted[col], rtol=1e-5, equal_nan=True):
                differences.append(col)
        else:
            # String comparison
            if not df1_sorted[col].equals(df2_sorted[col]):
                differences.append(col)
    
    if differences:
        print(f"  ❌ Differences in columns: {differences}")
        return False
    
    print("  ✅ Results match!")
    return True


def test_sequence(seq, name, linear=False):
    """Test a single sequence with both implementations."""
    print(f"\n{'='*60}")
    print(f"Testing: {name}")
    print(f"Sequence length: {len(seq)}")
    print(f"Linear: {linear}")
    
    try:
        # Run original implementation
        print("\nRunning original implementation...")
        original_results = original_annotate(seq, linear=linear)
        
        # Run new Snakemake implementation
        print("Running new Snakemake implementation...")
        with TemporaryDirectory() as tmpdir:
            new_results = new_annotate(
                seq, 
                linear=linear,
                output_dir=tmpdir,
                threads=4
            )
        
        # Compare results
        success = compare_dataframes(new_results, original_results)
        
        return success
        
    except Exception as e:
        print(f"  ❌ Error: {e}")
        return False


def main():
    """Run compatibility tests."""
    print("pLannotate Pipeline Compatibility Test")
    print("=" * 60)
    
    # Test sequences
    test_cases = [
        # (sequence, name, linear)
        (
            "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA",
            "Short linear CDS",
            True
        ),
        (
            "TCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTG",
            "pUC19 fragment (circular)",
            False
        ),
        (
            # A sequence with ori
            "TTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTT",
            "Origin fragment",
            False
        ),
    ]
    
    results = []
    for seq, name, linear in test_cases:
        success = test_sequence(seq, name, linear)
        results.append((name, success))
    
    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY:")
    print(f"{'='*60}")
    
    for name, success in results:
        status = "✅ PASSED" if success else "❌ FAILED"
        print(f"{name}: {status}")
    
    # Overall result
    all_passed = all(success for _, success in results)
    if all_passed:
        print("\n✅ All tests passed! The new pipeline is compatible.")
        return 0
    else:
        print("\n❌ Some tests failed. Please investigate.")
        return 1


if __name__ == "__main__":
    sys.exit(main())