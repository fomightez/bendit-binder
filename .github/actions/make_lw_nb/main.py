import os
with open("index.ipynb", 'r') as index_nb:
    lines = index_nb.readlines()
lines = "".join(lines)
lines = lines.replace('"lightweight_archive = False"','"lightweight_archive = True"')
with open("lw.ipynb", 'w') as output:
    output.write(lines)
