
import re
from pathlib import Path

file_path = Path(r'c:\Users\a-DJA\GIT\nat_zacros\nat_zacros\simulation_set.py')
content = file_path.read_text()

# Pattern: Header, then newline, then underline, then two or more newlines
pattern = re.compile(r'(Attributes|Parameters|Returns|Methods|Examples|Raises)\s*\n\s*(-+)\s*\n\s*\n')

for m in pattern.finditer(content):
    line_num = content.count('\n', 0, m.start()) + 1
    # Get some context
    lines = content.splitlines()
    context = "\n".join(lines[line_num-1 : line_num+5])
    print(f"Match found at line {line_num}:")
    print(context)
    print("-" * 20)
