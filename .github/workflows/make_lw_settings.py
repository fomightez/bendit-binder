import os
nb_to_change_settings_in = "index.ipynb"
with open(nb_to_change_settings_in, 'r') as index_nb:
    conent = index_nb.read()
content = conent.replace(
    '"lightweight_archive = False"','"lightweight_archive = True"')
with open(nb_to_change_settings_in, 'w') as output:
    output.write(content)
    
# adding extra code as extra cells based on 
# https://stackoverflow.com/a/45672031/8508004 and
# https://discourse.jupyter.org/t/delete-all-code-cells-except-markdown-text/3072/2?u=fomightez
import nbformat as nbf
ntbk = nbf.read(nb_to_change_settings_in, nbf.NO_CONVERT)
new_ntbk = ntbk

alertcode ="""\
import time
def speak(text):
    from IPython.display import Javascript as js, clear_output
    # Escape single quotes
    text = text.replace("'", r"\'")
    display(js(f'''
    if(window.speechSynthesis) {{
        var synth = window.speechSynthesis;
        synth.speak(new window.SpeechSynthesisUtterance('{text}'));
    }}
    '''))
    # Clear the JS so that the notebook doesn't speak again when reopened/refreshed
    clear_output(False)
speak('The run has completed, Nathan.')
for x in range(2):
    time.sleep(4)
    speak('The run has completed, Nathan.')"""

timecode = """\
import time
def executeSomething():
    #code here
    print ('.')
    time.sleep(480) #60 seconds times 8 minutes
while True:
    executeSomething()"""

new_ntbk.cells = ntbk.cells + [nbf.v4.new_code_cell(alertcode),
               nbf.v4.new_code_cell(timecode)]
nbf.write(new_ntbk, nb_to_change_settings_in, version=nbf.NO_CONVERT)
