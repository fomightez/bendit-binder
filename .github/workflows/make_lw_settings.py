import os
nb_to_change_settings_in = "index.ipynb"
with open(nb_to_change_settings_in, 'r') as index_nb:
    conent = index_nb.read()
content = conent.replace(
    '"lightweight_archive = False"','"lightweight_archive = True"')
with open(nb_to_change_settings_in, 'w') as output:
    output.write(content)
