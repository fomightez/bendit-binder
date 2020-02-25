import os
with open("index.ipynb", 'r') as index_nb:
    conent = index_nb.read()
content = conent.replace('"lightweight_archive = False"','"lightweight_archive = True"')
with open("lw.ipynb", 'w') as output:
    output.write(content)
