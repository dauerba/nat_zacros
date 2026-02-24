
import re
from pathlib import Path

file_path = Path(r'c:\Users\a-DJA\GIT\nat_zacros\nat_zacros\simulation_set.py')
lines = file_path.read_text().splitlines()

headers = {'Attributes', 'Parameters', 'Returns', 'Methods', 'Examples', 'Raises'}

for i, line in enumerate(lines):
    stripped = line.strip()
    if stripped in headers:
        # Check if next line is underline
        if i + 1 < len(lines):
            next_line = lines[i+1].strip()
            if next_line and all(c == '-' for c in next_line) and len(next_line) >= 3:
                # Found header and underline. Now check if next line is empty.
                if i + 2 < len(lines):
                    after_underline = lines[i+2].strip()
                    if not after_underline:
                        print(f"Match found at line {i+1}:")
                        print(f"{i+1}: {lines[i]}")
                        print(f"{i+2}: {lines[i+1]}")
                        print(f"{i+3}: {lines[i+2]}")
                        if i + 3 < len(lines):
                            print(f"{i+4}: {lines[i+3]}")
                        print("-" * 20)
