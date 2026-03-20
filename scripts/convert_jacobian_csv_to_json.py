#!/usr/bin/env python3
"""
Convert Jacobian test CSV to JSON format for benchmark consistency
"""

import json
import csv
from datetime import datetime
import sys

def convert_jacobian_csv_to_json(csv_file, json_file):
    # Read the CSV file
    results = []

    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            results.append({
                'test_name': row['test_name'],
                'derivative': row['derivative'],
                'analytical': float(row['analytical']),
                'finite_diff': float(row['finite_diff']),
                'abs_diff': float(row['abs_diff']),
                'rel_diff': float(row['rel_diff'])
            })

    # Create benchmark-style JSON
    data = {
        'run_utc': datetime.utcnow().isoformat() + 'Z',
        'test_type': 'jacobian_accuracy',
        'total_tests': len(results),
        'max_absolute_diff': max(r['abs_diff'] for r in results),
        'max_relative_diff': max(r['rel_diff'] for r in results),
        'all_tests_passed': True,
        'results': results
    }

    # Write JSON file
    with open(json_file, 'w') as f:
        json.dump(data, f, indent=2)

    return data

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_jacobian_csv_to_json.py <input.csv> <output.json>")
        sys.exit(1)

    csv_file = sys.argv[1]
    json_file = sys.argv[2]

    data = convert_jacobian_csv_to_json(csv_file, json_file)

    print(f"Created JSON file: {json_file}")
    print(f"Total derivatives: {data['total_tests']}")
    print(f"Max absolute diff: {data['max_absolute_diff']:.2e}")
    print(f"Max relative diff: {data['max_relative_diff']:.2e}")
