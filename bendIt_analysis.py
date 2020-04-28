# bendIt_analysis.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"


# bendIt_analysis.py by Wayne Decatur
# ver 0.1.0
#
#*******************************************************************************
# 
# PURPOSE: Processes via standalone bendIt many 'cassette' sequences placed in 
# defined flanking sequences.
# Needs to work in conjunction with the notebook `index.ipynb` that is present
# at https://github.com/fomightez/bendit-binder . In fact, the easiest way to
# use this is to launch sessions by clicking on the `launch bend.it` badge at
# that repo. In the session that comes up, everything will already be installed
# and presented to the user for processing.
# Originally written as an Ipython script, `bendIt_analysis.ipy` v.0.2.2, it has 
# been converted to pure Python. Offers considerable speed advantage too, 
# (50 minutes originally down to 9.5 for 658 sequences of around 75 bp total) 
# probably because not making all those separate shells for file deletions, etc.
# 
import os
import sys
import glob
from shutil import copyfile
import subprocess
import pandas as pd
import pyfaidx
from halo import HaloNotebook as Halo
from IPython.utils import io
import fnmatch
import datetime
from statistics import mean, median, mode, StatisticsError
import collections

demo_file_name = "demo_sample_set.fa"

fasta_extensions = ("fa", "fasta", "faa", "fas", "fsa")



if lightweight_with_images:
    # use lightweight settings but also make plot images
    lightweight_archive = True
    make_images = True
elif not lightweight_archive:
    # use standard settings that includes making plot images
    make_images = True
else:
    make_images = False
    



################################################################################
#######----------------------HELPER FUNCTIONS-----------------------------######

def write_string_to_file(s, fn):
    '''
    Takes a string, `s`, and a name for a file & writes the string to the file.
    '''
    with open(fn, 'w') as output_file:
        output_file.write(s)

def out2_stderr_n_log(s,log_file_text):
    '''
    Takes a string as input and sends it to the stderr as well as to a building
    string that will everntually get saved as a Log file.
    Also needs the Log file to be sent in because gets assigned within the
    function in order to add to it. Returns the modified `log_file_text`.
    '''
    sys.stderr.write(s)
    log_file_text += s
    return log_file_text


import time
from IPython.display import display, Javascript
import hashlib
def save_notebook(file_path):
    '''
    Function to save a notebook from 
    https://stackoverflow.com/a/57814673/8508004

    IMPORTANTLY, this won't work in the JupyterLab interface for notebooks!
    See https://github.com/jupyterlab/jupyterlab/issues/7627
    '''
    start_md5 = hashlib.md5(open(file_path,'rb').read()).hexdigest()
    display(Javascript('IPython.notebook.save_checkpoint();'))
    current_md5 = start_md5
    while start_md5 == current_md5:
        time.sleep(1)
        current_md5 = hashlib.md5(open(file_path,'rb').read()).hexdigest()


def sanitize_description_lines(fn,delimiter_character):
    '''
    User wants to be able to have slashes for part of a date stamp in the 
    description lines. This scans FASTA file for description lines and removes
    those slashes as they will interfere with using those names as tags in file
    names. (example description line with problem: `R5|11/23/19|16(4) `)
    While at it going to remove and trailing whitspace at end so I can easily 
    use the sample name derived from description line for tracking and not have
    to worry about removing space for ease in use in file names and lose 
    correspondence to actual pyfaidx keys for each sequence.

    Added in also taking in `delimiter_character` and using that to limit what
    needs sanitizing. Don't care if it is outside of what will become sample 
    name.

    Extra sanitizing added when I found bendIt/bash didn't work with `|` or 
    parantheses in the file names which I was using the description line as.
    Will use position of closing paranthesis at end to get away using same 
    indicator character for both since running out of unusual characters to use.
    Use of `+` from https://superuser.com/questions/358855/what-characters-are-safe-in-cross-platform-file-names-for-linux-windows-and-os
    (I wonder if these listed at 
    https://twitter.com/PhilippBayer/status/1229741303727935493, "ASCII has a 
    set of characters especially for this purpose. Nobody uses it.  ASCII 1C to 
    1F" would be useable?!?! It would add more options that I had come up with 
    for replacements that I could the substitute back out later.)
    Example shown above gets sanitized to `R5_11-23-19_16+4+`

    Takes:
    a file name

    Returns:
    Nothing but as part of process saves clean version of file as same name as
    file name argument. NOTE: If nothing is encountered that needs cleaning, 
    nothing is done to input file.
    '''
    anything_needing_sanitizing_encountered = False
    collecting_content_string = ""
    with open(fn, 'r') as input:
        for line in input:
            if line.startswith(">"):
                line = line.split(delimiter_character)[0]
                if "/" in line:
                    line = line.replace("/","-")
                    line = line.replace("|","_")
                    line = line.replace("(","+")
                    line = line.replace(")","+")
                    anything_needing_sanitizing_encountered = True
                    line = line.strip() + "\n" #might as well remove any trailing 
                    # space as well while handling slashes because no harm if 
                    # moot & then don't need to check that line again in any way 
                    # to rule if possesses trailing space. ACtually when added
                    # in `delimiter_character` to set as boundary to sample name
                    # the `strip` became MOOT when default space is used.
                    collecting_content_string += line
                elif line.endswith(" "):
                    line = line.strip() + "\n" #remove trailing space
                    #anything_needing_sanitizing_encountered = True #Changed 
                    # from this being noted as 
                    # `anything_needing_sanitizing_encountered` when added 
                    # `delimiter_character` use because default is space. And
                    # if using that default end in space will be ignored anyway,
                    # and if using someting else should hit before a space at 
                    # end.
                    line = line.strip() + "\n" #remove trailing space
                    collecting_content_string += line
                else:
                    collecting_content_string += line
            else:
                    collecting_content_string += line
    # replace input file if any sanitizing occured.
    if anything_needing_sanitizing_encountered:
        #!cp fn unsanitized_{fn} #cannot use in function, like `%store` because
        # actually run in global namespace which cannot see Python variable `fn`
        # local to the function
        sanitized_fn_name = f"unsanitized_input_{fn}.txt"
        copyfile(fn, sanitized_fn_name) #adding extension means it won't get 
        # recognized as a FASTA if script is run again
        sys.stderr.write("The original input was saved as '{}' because "
            "characters\nin description line needed replacing in order to\nbe "
            "properly processed by this script & bendIt.".format(
            sanitized_fn_name))
        write_string_to_file(collecting_content_string, fn)

def chunk_string(string, chunk_size):
    """Return a list of n-sized chunks from string of letters."""
    return [string[i:i+chunk_size] for i in range(0, len(string),chunk_size)] 


def strip_off_first_line(fn,set_name,character_to_mark_set_name_end):
    '''
    This takes a name of a file & then uses the shell to remove the first line.
    In order to leave the input file intact, a new multi-sequence FASTA file
    is made and that is used in place of the one where the label was the first
    line. The set sample name extracted gets added to the file name.
    Removing first line based on 
    https://unix.stackexchange.com/questions/96226/delete-first-line-of-a-file
    '''
    name_for_f_without_first_line = (
        f"{set_name}{character_to_mark_set_name_end}set.fa")
    #!tail -n +2 {fn} >{name_for_f_without_first_line} 
    os.system(f"tail -n +2 {fn} >{name_for_f_without_first_line}")
    return name_for_f_without_first_line

def divide_up_any_additional_sets(f,character_to_mark_set_name_end):
    '''
    Takes a file that may contain additional sample set labels and fasta entries
    that correspond and makes files for each. It returns the list of files made.

    This only works on a file where the first label in the first line has
    already been removed so the first line in the current file starts with a 
    description line. It relies on the sequences in the input file requested
    by the user being short (or at least all on one line) and the fact the
    inclusion of the labelwould disrupt the cycle that starts with the
    description line and then is a sequence line and then another description
    line and sequence. The provided cassette sequences should continue on that
    way unless a label for a sample set occurs. In that case, a new set is
    started. If another label is then later encountered the current string is
    saved to a file and that cycle of looking for another label to designate
    another set repeats on again.
    '''
    # first copy the input file that will be parsed line by to a new file so
    # can parse contents while possibly overwriting the input file with a
    # shorter version if a label for a set encountered inside it
    temp_file_name = "temp_copy.fa"
    #!cp {f} {temp_file_name}
    copyfile(f, temp_file_name)
    # set up some variables for holding assignments as go through line by line
    additional_sequence_files = []
    current_string = ""
    current_set_name = "MOOT" #moot at this point b/c already in name
    current_file_to_write_to = f
    # go through the file line by line until cycle of description line following 
    # sequence line is disrupted
    with open(temp_file_name, 'r') as input:
        nxt_line_expct_descript, nxt_line_expct_seq = (True,False)
        for line in input:
            if nxt_line_expct_descript and (
                line.startswith(">") == False) and line.strip() != "":
                # label encountered; so write the current string that was being
                # built and set variables for collecting contents going forward
                write_string_to_file(current_string, current_file_to_write_to)
                # If there were internal labels, need to record the just saved
                # file to the additional file list, but not if it was the
                # original file that was provided when the function called.
                if current_file_to_write_to != f:
                    additional_sequence_files.append(current_file_to_write_to)
                current_set_name = line.strip()
                current_file_to_write_to = (
                    f"{current_set_name}{character_to_mark_set_name_end}set.fa")
                current_string = ""
                nxt_line_expct_descript, nxt_line_expct_seq = (True,False)
            elif nxt_line_expct_descript and line.startswith(">"):
                current_string += line
                nxt_line_expct_descript, nxt_line_expct_seq = (False,True)
            elif nxt_line_expct_seq and (line.startswith(">") == False):
                current_string += line
                nxt_line_expct_descript, nxt_line_expct_seq = (True,False)
            elif line.strip() == "":
                #there is no need to do anything if the line is blank like could 
                # be at end of document or possibly(?) sample set section. The 
                # user supplied examples have not shown evidence of the latter.
                continue
            else:
                issue = (
                    "This shouldn't happen since conditions checked should "
                    "catch all possibilies.")
                sys.stderr.write(f"\n{issue}")
    # If no internal labels encountered then the original file won't need to be
    # divided (shortened) at all and so there is no need write anything. On the
    # other hand, if there were internal labels, this next line will save
    # the current string being built to a fasta file for the current set. And if
    # there were internal labels, need to record just saved file to the
    # additional file list.
    if current_file_to_write_to != f:
        write_string_to_file(current_string, current_file_to_write_to)
        additional_sequence_files.append(current_file_to_write_to)
    #clean up by deleting the copy of the input that was used to iterate line by
    #line while possibly modifying original fasta file.
    #!rm {temp_file_name}
    os.remove(temp_file_name)
    return additional_sequence_files

def percent_GCcalc(items):
    '''
    takes a list of three and calculates percentage of sum of first
    two itemswithin total (second item)

    Taken from 
    `GSD Adding_percentGC_to_nt_counts_for_mito_genomes_from_1011_collection.ipynb`
    '''
    return (items[0] + items[1])/items[2]

def make_and_run_review_nb(now, review_nb_stub, serial_fn):
    '''
    Takes current datetime, a string that represents the notebook 'stub' (in a 
    python script form that later is converted to a proper nb by Jupytext), and 
    the file name of the serialized data and makesa Python script. It then 
    converts that script to a notebook and runs it using jupytext. The produced 
    notebook contains every plot displayed for review.

    It returns the name of the notebook produced
    '''
    plots4review_fn = f"plots4review_from-ba{now.strftime('%b%d%Y%H%M')}.ipynb"
    review_nb_stub = review_nb_stub.replace(
        "PICKLED_SEQS_DFS_N_PLOTS_FILE_PLACEHOLDER",serial_fn)
    write_string_to_file(review_nb_stub, plots4review_fn[:-6]+".py")
    #!jupytext --to notebook --execute {plots4review_fn[:-6]+".py"}
    os.system(f'jupytext --to notebook --execute {plots4review_fn[:-6]+".py"}')
    return plots4review_fn

#######------------------END OF HELPER FUNCTIONS--------------------------######
################################################################################


# delete `start` to make the starting directory cleaner since will be messy
# while running & had wanted to delete it as part of `start` but didn't want to
# bother seeing if possible.
for file in os.listdir('.'):
    if file == "start":
        #!rm start
        os.remove("start")



#For use with converting `report_with_curvature` settings back to what
# bendIt needs

'''
from enum import Enum
class report_with_curvature_correspondences(Enum):
    G+Ccontent = "G"  #<--- Cannot use Enums because I think the plus will mess things up
    bendability = "B"
    complexity = "C"
'''
report_with_curvature_settings_corrspndnce = {
    "G+Ccontent" : "G",
    "bendability" : "B",
    "complexity" : "C"
}

#So in command it would be like:
# !bendIt -s test.fa -o result.txt -g {report_with_curvature_settings_corrspndnce[report_with_curvature]}

# a similar correspondance for giving dataframe columns names
results_to_df_corrspndnce= {
    "G+Ccontent" : "GCcontent",
    "bendability" : "Bendability",
    "complexity" : "Conplexity"
}

curvature_window_size,other_metric_reporting_window_size =(
window_size,window_size) #because the window size gets used for short sequences 
# to determine what to plot, it is best if these are linked and the same. In 
# development, I was letting them be set different
# early on like the bend.it server does but this script does some handling 
# behind-the-scenes to get things to work with short sequences, and best to lock
# into single value, and so this insures now they are, but doesn't make me 
# change code used below here.

log_file_text = ("####BELOW IS LOG of std.out and std.err from bendit "
    "analysis: ######\n\n")


review_nb_stub ='''# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Plots resulting from bendIt Analysis for review
#
# For ease of review, use JupyterLab to place 'views' of this notebook side-by-side. As many as seem reasonable to fill your screen width.

import pickle
serial_fn = "PICKLED_SEQS_DFS_N_PLOTS_FILE_PLACEHOLDER"
with open(serial_fn, "rb") as f:
    seqs_dfs_and_plots_per_sample_set = pickle.load(f)


# %matplotlib inline
def make_manager(fig):
    # create a dummy figure and use its
    # manager to display "fig"  ; based on https://stackoverflow.com/a/54579616/8508004
    dummy = plt.figure()
    new_manager = dummy.canvas.manager
    new_manager.canvas.figure = fig
    fig.set_canvas(new_manager.canvas)
import seaborn as sns
import matplotlib.pyplot as plt
for ss,collected_dicts in seqs_dfs_and_plots_per_sample_set.items():
    for sample_name,the_plot in collected_dicts[-1].items():
        make_manager(the_plot)
        the_plot.show()

''' 









################################################################################
#######------------------------MAIN SECTION-------------------------------######

spinner = Halo(text='Processing...', spinner='dots',color = 'magenta')
spinner.start()


# GO THROUGH AND SET-UP SEQUENCE FILES TO ANALYZE:
#------------------------------------------------------------------------------#
# This section prepares the sequence files to analyze if they are non-standard,
# i.e., have labels inside what is called the FASTA file. An input option 
# requested by a user
# This section also determines if the demo sequences are to be used and does 
# some checks even if input file(s) are standard FASTA. Part of that is seeing
# if extracting the sample set names or determining from file name.
sequence_files = []
for file in os.listdir('.'):
    if file.endswith(fasta_extensions):
        sequence_files.append(file)
# Disregard the demo file if any other sequence files added
if len(sequence_files) > 1 and demo_file_name in sequence_files:
    sequence_files.remove(demo_file_name)
#Report if demo is being used or how many sequence files being processed
if len(sequence_files) == 1 and sequence_files[0] == demo_file_name:
    for_out2both = "Processing demo file '{}'.\n".format(demo_file_name)
    log_file_text = out2_stderr_n_log(for_out2both,log_file_text)
else:
    for_out2both = ("Processing {} files: {}.\n".format(
        len(sequence_files),repr(sequence_files)))
    log_file_text = out2_stderr_n_log(for_out2both,log_file_text)

if lightweight_archive:
    lw_note_text = ("\n'lightweight_archive' setting used. The resulting "
        " archive will not\n include any 'raw' gnuplots, notebook for review "
        "of plots, or image files.\n The latter two can be made later from the "
        "archive, if needed. The LOG file\nwill not report on individual "
        "sequence processing.\n")
    log_file_text = out2_stderr_n_log(lw_note_text,log_file_text)

# Go through list of sequence_files and check if these are standard FASTA or if
# combine sample set label(s) and FASTA record(s). Although the most noticeable 
# feature of this situation is that there is a label at the first line of the
# input file followed by a description line, the label for an additional sample
# set can occur later as well, and so this then requires extra special handling.
# Namely each sample set has to be saved to a standard multi-seq FASTA file, 
# with any internal sets requiring still another standard FASTA file made 
# corresponding to that set. This handling is going to rely on enforcing that 
# for these non-standard, multi-seq FASTA files, the sequences for the cassettes
# must be constrained to a single line, or it will be too hard to discern in the 
# inernal cases if a next line is a label or sequence. The list of 
# `sequence_files` also needs adjusting upon any special handling to standardize 
# the input to contain the updated names of the files where the label is in the
# file name for the set.
label_in_fasta = False
file_name_substitutions_to_make = {} # this will be used if non-standard input 
# file with label on first line is used
additional_sequence_files = [] # This will get used if non-standard input file #
# used with labels more interal beyond just the first line
for x in sequence_files:
    # Sanitize description lines to accommodate user wanting to include date 
    # stamps with forward slashes as part of FASTA description lines. 
    # Sanitizing will make it possible to use this as part of a file name 
    # more easily in order to easily track sequences and results throughout.
    # Also going to remove line-ending spaces while at it just to be cleaner.
    sanitize_description_lines(x,character_to_mark_individual_sample_name_end)
    # determine if label for the name of the sample set is on first line or if
    # using file name or defined `sample_set_name`
    if sample_set_name_extract_auto:
        #check if label for the name of the sample set is on first line by
        # checking if first `>` on second line
        #first_two_lines_list = !head -2 {x} #when use that it goes to a string 
        # and not bytes but the subprocess way produces byets so changed to add 
        # conversion to string so didn't need to change any other code below
        first_two_lines_list = subprocess.check_output(
            f"head -2 {x}", shell=True).split()
        first_two_lines_list = (
            [ft.decode("utf-8") for ft in first_two_lines_list]) #the subprocess
        # way produces byets so changed to add conversion to string so didn't 
        # need to change any other code below
        if ((not first_two_lines_list[0].startswith('>')) and (
            first_two_lines_list[1].startswith('>'))):
            label_in_fasta = True
            sample_set_name = first_two_lines_list[0].strip()
            without_first_line_f = strip_off_first_line(
                x,sample_set_name,character_to_mark_set_name_end)
            file_name_substitutions_to_make[x] = without_first_line_f
            # Because if the sequence file has been determined to be 
            # non-standard & have labels in the file, there is the possibility
            # additional sample set labels exist later in the file beyond the 
            # first line & additional sections should be split out into 
            # separate, standardized fasta files.
            additional_sequence_files = divide_up_any_additional_sets(
                without_first_line_f, character_to_mark_set_name_end)
        else:
            sample_set_name = x.split(character_to_mark_set_name_end)[0].strip()
    # 'else' situation will be that `sample_set_name` is used
    else:
        assert len(sequence_files) == 1, "To use "
        "`sample_set_name_extract_auto = False` there can only be one input "
        "file as keeping the correspondences intact would be cumbersome for "
        "multiple input files and multiple name settings."
        assert sample_set_name, "`sample_set_name` must be assigned as the name "
        "to use for the sample set if `sample_set_name_extract_auto` is set to "
        "`False`."
        # Later in script, handle using the user-defined `sample_set_name` for
        # the results. Because other obvious choice would be to add name into 
        # FASTA file name and I'd have to be sure it was compatible with 
        # getting it back out. Might as well just make it easier.
# If non-standard input FASTA files were used, then the list of `sequence_files`
# may need adjusting.
if label_in_fasta:
    for old,new in file_name_substitutions_to_make.items():
        sequence_files = [f.replace(old,new) for f in sequence_files]
    sequence_files += additional_sequence_files

# GO THROUGH SEQUENCE FILES TO ANALYZE:
#------------------------------------------------------------------------------#
# Now step through using the standardized FASTA files. In most cases, the sample 
# names will be built in to the file name; however, there is the possibility 
# that for single input files, the option  of setting 
# `sample_set_name_extract_auto = False` and assigning the name with 
# `sample_set_name`.
#files_produced = []
files_produced = sequence_files.copy() #start by including sequence files so 
# they are stored along with that output in a single archive at the end
seqs_dfs_and_plots_per_sample_set = {} #overarching collection keyed by sample 
# set; each set of items keyed by sample will go into this as a list
dataframes_produced_keyed_by_sample = {}
plots_produced_keyed_by_sample = {}
for indxf,x in enumerate(sequence_files):
    # determine if label for the name of the sample set is on first line or if
    # using file name or defined `sample_set_name`
    label_in_fasta = False
    if sample_set_name_extract_auto:
        #check if label for the name of the sample set is on first line by
        # checking if first `>` on second line
        #first_two_lines_list = !head -2 {x}
        first_two_lines_list = subprocess.check_output(
            f"head -2 {x}", shell=True).split()
        first_two_lines_list = (
            [ft.decode("utf-8") for ft in first_two_lines_list]) #the subprocess
        # way produces byets so changed to add conversion to string so didn't 
        # need to change any other code below
        if ((not first_two_lines_list[0].startswith('>')) and (
            first_two_lines_list[1].startswith('>'))):
            label_in_fasta = True
            sample_set_name = first_two_lines_list[0].strip()
            x = strip_off_first_line(x)
        else:
            if character_to_mark_set_name_end in x:
                sample_set_name = x.split(
                    character_to_mark_set_name_end)[0].strip()
            # Build in some fall backs for extracting name of the set, in case 
            # advice about naming isn't followed.
            elif " " in x:
                sample_set_name = x.split(" ")[0]
            elif "-" in x:
                sample_set_name = x.split("-")[0].strip()
            else:
                sample_set_name = x.split(".")[0].strip()

    # 'else' situation will be that `sample_set_name` is used
    else:
        assert sample_set_name, "`sample_set_name` must be set as the name to "
        "use for the sample set if `sample_set_name_extract_auto` is set to "
        "`False`."
    # Iterate over each sequence to place in the cassette position with defined
    # flanking sequence and use bendIt to analyze . For ease of reference, the
    # sequences withe the cassette places in the flanking sequences will be
    # called 'merged' here.
    sequences_processed_keyed_on_name = {} # also for keeping track of 
    # subsequently without need for going back to FASTA
    merged_sequences_processed_keyd_on_name = {} # for tracking the pieced
    # together final sequences
    sequences_with_duplicates = False
    merged_sequences_with_duplicates = False
    from pyfaidx import Fasta
    sequence_entries = Fasta(x, read_long_names=True)
    #print(len(sequence_entries))  # FOR DEBUGING ONLY
    for indxe,sequence in enumerate(sequence_entries):
        sample_name = sequence.name.split(
            character_to_mark_individual_sample_name_end)[0].strip()
        sequences_processed_keyed_on_name[sample_name] = str(sequence).strip()
        merged_sequence = f"{up_bracket}{str(sequence).strip()}{down_bracket}"
        merged_sequences_processed_keyd_on_name[sample_name] = merged_sequence
        # make a fasta file (with a single sequence) out of the merged sequence
        fasta_file_name_for_merged = f"{sample_name}_merged.fasta"
        # If merged sequence is actually less than 1019 (although found in other 
        # testing the exact number varies and sometimes will be around 1040 or 
        # 1030) then add repeats of the sequence to make it longer than that or 
        # otherwise get 
        # 'Segmentation fault (core dumped)' because it seeme the value at the
        # 5th base goes to infinity and crashes the process. Doesn't happen with
        # longer sequences and if clean results back to the length of the window
        # it shouldn't matter. Going to round up a bit from 1019 since it seemed
        # the exact cutoff differed in different sessions as it could possibly
        # be tied to available computing power what cutoff raises the issue.
        # Was originally using dummy sequence of a string of `T`s to add in turn
        # to each end but was finding result for predicted curvature value for
        # last position that gives a curvature value other than zero wasn't 
        # matching what bend.it server gives for that position. Also by testing
        # addition of strings of different nucleotides of same length was noting
        # the identity the nucleotides of the string added affected result even
        # though number of nucleotides added same. So sequence somehow effects
        # it and upon testing adding repeats I found for sequences of varying
        # length that I could add certain multiples of repeats and then get the
        # value for the last curvature reporting position the same as bend.it
        # server reports. This was actually worked out with sequences long 
        # enough to not throw segmentation error so I could track a series of 
        # repeats, starting with none (just the input sequence) up to around 50
        # as a I lengthened the sequence. Interestingly depedning on length of
        # the sequence I would see different multiples of repeats work. For 
        # example, one sequence always worked but I found the addition of 21,7,
        # or 3 repeats seemed to work depending on length. Interstingly, all
        # these intersect at 21 and multiples of 21 (42 is only one I check but 
        # assume beyond) work. So that is probably easiest to work with instead
        # of determinr best value for each length below the cutofff. I believe
        # the 21 magic number comes from there being 10.5 bases per helical turn
        # and so the effect of the previous sequence on the curvature of the 
        # last position reporting curvature is the same when it is phased into
        # the same relative position of the sequence with the same helical
        # surface properties proceeding it. Note that I did sometimes see the
        # curvature value for this last position vary by 0.0001 from what the
        # bend.it server reported for 21 repeats but never more. This is a 
        # respectable deviation from what the bend.ot server reports because I 
        # can say the value for that last position agrees always to a thousandth 
        # of a decimal and actually never varies by more than 0.0001. And that 
        # seemsreasonable deviaiton to accept for convenience of being able to 
        # implement this with one repeat value.

        # Because need to know the size of merged_sequence without line breaks,
        # and also before any additional sequence POSSIBLY added, going to
        # to first assign the string that will get linebreaks direct from 
        # `merged_sequence` and then I can possibly add to that string if it is
        # still too short without altering the true `merged_sequence`
        merged_sequence_with_linebreaks = merged_sequence
        cutoff_for_adding = 1041
        if len(merged_sequence_with_linebreaks) < cutoff_for_adding:
            # if the sequence alone is not long enough to avoid a segmentation
            # error when bendIt command is run, add repeats of the sequence to
            # insure it is long enough and that the far right end can be mined
            # for what the sequence at the end would give if no additional 
            # sequence added
            while len(merged_sequence_with_linebreaks) < cutoff_for_adding:
                merged_sequence_with_linebreaks += merged_sequence*21
            sequence_repeated_to_not_cause_fault = True
        else:
            sequence_repeated_to_not_cause_fault = False
        merged_sequence_with_linebreaks = "\n".join(
            chunk_string(merged_sequence_with_linebreaks,70)) # `t.fasta` that 
        # came in the bendit `test` folder has 70 bases per line and in some 
        # early tests it seemed if I didn't have 70 bp per line I got 
        # 'Segmentation fault (core dumped)'; however, then I check really long 
        # sequences later it seemed it didn't matter ACTUALLY. Left chunking to 
        # 70 bp in just because it makes it easier to tell length in developing.
        merged_sequence_made_to_fasta = (
            f">{sample_name}\n{merged_sequence_with_linebreaks}")
        # keeping this commented out for now because here `%store is working
        #!echo {merged_sequence_made_to_fasta} > {fasta_file_name_for_merged} 
        # above uses shell b/c magics like `%store` work poor when in a function 
        # and considering may move to a function later
        #Switching from using `%store` command seen below to this line when bash 
        # was choking on file names like `R5|11-1-19|1(1)_merged.fasta` that was 
        # derived from what a user needed in the description line to indicate 
        # the sample.
        write_string_to_file(
            merged_sequence_made_to_fasta, fasta_file_name_for_merged)
        #%store merged_sequence_made_to_fasta >{fasta_file_name_for_merged} #it 
        # is important to not have a space after  redirect; apparently recent 
        # change (because I only noted in Jan. 2020) now reflected in 
        # documentation at 
        # https://ipython.readthedocs.io/en/stable/config/extensions/storemagic.html
        
        # Use the merged sequence and the settings to run bendIt.
        # Note that in order to get the standalone version to work and not give
        # give a segmentation fault and dump the core, I had to use multiples
        # of the input sequence in cases it was shorter than around 1040. The 
        # bendIt command always produces a gnuplot and to make it match mostly
        # to the input merged squence can use  --xmin` & --xmax` options so
        # plot is at least for range I care about. I don't plan to use these 
        # plots for anything but review in development but seems could be useful
        # for development for them to plot the data of interest.
        end_range = len(merged_sequence) - curvature_window_size
        output_file_suffix = f"{sample_name}_output"
        with io.capture_output() as captured:
            '''
            !bendIt -s {fasta_file_name_for_merged} -o {output_file_suffix} -c \
            {curvature_window_size} -b {other_metric_reporting_window_size} \
            -g {report_with_curvature_settings_corrspndnce[report_with_curvature]} \
            --xmin {curvature_window_size} --xmax {end_range}
            '''
            os.system(f"bendIt -s {fasta_file_name_for_merged} -o \
                {output_file_suffix} -c {curvature_window_size} -b \
                {other_metric_reporting_window_size} -g \
                {report_with_curvature_settings_corrspndnce[report_with_curvature]} \
                --xmin {curvature_window_size} --xmax {end_range}")

            # POST PROCESS RESULTS:
            #------------------------------------------------------------------#
            # Bring the reported values into Python for convenient plotting & 
            # analyzing & storing.

            #check for script to send bendit standalone results to dataframe(s)
            # and retrieve if not already here
            # Get a file if not yet retrieved / check if file exists
            file_needed = "bendit_standalone_results_to_df.py"
            if not os.path.isfile(file_needed):
                #!curl -OL https://raw.githubusercontent.com/fomightez/sequencework/master/bendit_standalone-utilities/{file_needed}
                os.system("curl -OL https://raw.githubusercontent.com/"\
                    "fomightez/sequencework/master/"\
                    f"bendit_standalone-utilities/{file_needed}")

            # convert to a dataframe
            from bendit_standalone_results_to_df import bendit_standalone_results_to_df
            df = bendit_standalone_results_to_df(
                output_file_suffix,
                results_to_df_corrspndnce[report_with_curvature],
                sequence_file=fasta_file_name_for_merged,
                pickle_df=False) # Don't pickle in external script call because 
            # want to pickle and save tab-delimited files here where can incoporate 
            # the sample name & add the dataframe and and saved files to a 
            # collection of files produced. Plus, what is made when short sequences
            # is artifically long to get past segmetnation fault issue & neds to be 
            # fixed at right-had end and so that is a further reason to wait and 
            # handle anything along those lines here.
        log_file_text += captured.stdout
        log_file_text += captured.stderr

        # If the sequence had to be repeated multiple times because it otherwise
        # would have been short and caused asegmentation fault when bendIt 
        # runs, the right-hand side of the initial sequence will be affected by 
        # the addiitonal repeats downstream. To get data for that end that isn't 
        # affected by additional sequence, the far end on the right, i.e., the 
        # end of the last repeat will be used to substitute in for that end. And 
        # because the right number of repeats was used the value for the last 
        # reported curvature will be as it would be if only the inut sequence 
        # has been used because the effect from the side of the helix will
        # match. Note this ends up also being more efficent then when the dummy 
        # sequence was being used because in tht case I used to add sequence
        # at  the left-hand side (upstream end/ 5'-end) to avoid segmentation 
        # fault and run that through the bendIt command again and process the 
        # datat from that round to get data I thought would be useful for the 
        # right-hand side. (Turns out it was good for most values but not for 
        # the last reported curvature value.)
        # This route to avoid the segmentation fault has been verified to give 
        # the same results as the bend.it server to within a thousandth of a 
        # decimal with numerous sequences. In fact, I never saw more than a 
        # 0.0001 difference, with any difference being a very rare occurence.
        # Substituing in the data from the far right hand side to the end of 
        # the input seqeunce.
        if sequence_repeated_to_not_cause_fault:
            # Set up merging it to the left-hand side of original input sequence
            point_to_merge_at = (end_range - ((curvature_window_size + 3) + 2))#
            # The `+ 3` portion comes from the fact bendability is calculated for 
            # trinucloetide values as dicussed at http://pongor.itk.ppke.hu/dna/bend_it.html#/bendit_theory 
            # & shown in tables at http://pongor.itk.ppke.hu/dna/bend_it.html#/bendit_tables
            # and so on tip of the window size issue cutting off all values 
            # below the window size, the bendability isn't also calculated until
            # another trinucelotide has occured beyond the window size into the
            # sequence at the start and also at the far end. The latter end 
            # being the pertinent one here. It is '+ 2' beyond b/c going 
            # just beyond the first trinucleotide by one base didn't seem to get 
            # away from end effect entirely for some reason but a one more did 
            # so when at that point the curvature values in the dataframe from 
            # the first round and where sequence added at other end now same.
            #print("point_to_merge_at:",point_to_merge_at)  # FOR DEBUGING ONLY
            point_to_merge_at_far_end = (len(df) - (
                (curvature_window_size + 3) )) # don't need to subtract 
            # window from length of df because the df comes from bendIt results 
            # which have removed it.  I think I don't also need the +2 because
            # of zero indexing relative the length and want one after the merge 
            # point for the right side.
            # Verify the values at the merge point are indeed the same. They 
            # should not diverge until slightly farther down.
            assert (
                df.iloc[point_to_merge_at].Predicted_curvature == 
                df.iloc[point_to_merge_at_far_end].Predicted_curvature),("The "
                "values for curvature at the merge point should be the same.")
            assert (
                df.iloc[point_to_merge_at][df.columns[-1]] == 
                df.iloc[point_to_merge_at_far_end][df.columns[-1]]),("The "
                "values at the merge point should be the same.")
            # Merge the parts
            df=pd.concat(
                [df[:point_to_merge_at],df[point_to_merge_at_far_end:]], axis=0)
            df = df.reset_index(drop=True) # rest the index because don't need
            # the high index values after the mege point showing up and the 
            # fixed index values can be used to fix the positions
            # Need to fix the position numbers after the merge point. Can just 
            # make all set relative the index and it will fix the end ones too.
            df.Position = df.index + curvature_window_size





        # Plot the reported data
        #----------------------------------------------------------------------#
        spinner.stop()  #turn off spinner so it doesn't make only the most 
        # recent plot show
        import seaborn as sns
        import matplotlib.pyplot as plt
        sns.set_style("whitegrid")
        fig, ax1 = plt.subplots(figsize=(8,4))

        if report_with_curvature == "bendability":
            #ax2 = ax1.twinx() #Nathan is plotting bendability and curvature on same y-axis which won't work with the other two since their range is so low
            '''
            sns.lineplot(x=df.Position,
                        y=df.Predicted_curvature,
                        color="C0",
                        ax=ax1)
            sns.lineplot(x=df.Position, 
                         y=df.Bendability,
                         color="C1",
                         ax=ax2)
                         '''
            if smooth_plot_curves:
                # use of Akima1DInterpolator based on UPDATE#2 at 
                # https://stackoverflow.com/a/37223558/8508004 and seems to match
                # fairly well to Excel's 'Scatter with Smooth Lines' option
                import numpy as np
                from scipy.interpolate import Akima1DInterpolator
                akima1 = Akima1DInterpolator(df.Position, df.Predicted_curvature)
                akima2 = Akima1DInterpolator(df.Position, df.Bendability)
                x_4smooth = np.linspace(
                    min(df.Position), max(df.Position), round(
                    len(df)*14.492753)) # `len(df)*14.492753)` is about 1000 
                # when df is from 75 bp input (69 total accounting for 
                # window size)
                sns.lineplot(x=x_4smooth,
                            y=akima1(x_4smooth),
                            label="Curvature",
                            linewidth=1.9,
                            #color="C0",
                            color="xkcd:cerulean", #see 
                            # https://matplotlib.org/tutorials/colors/colors.html & 
                            # https://xkcd.com/color/rgb/
                            ax=ax1)
                sns.lineplot(x=x_4smooth, 
                             y=akima2(x_4smooth),
                             label="Bendability",
                             linewidth=1.9,
                             #color="C1",
                             color="xkcd:orange", #see 
                            # https://matplotlib.org/tutorials/colors/colors.html &
                            # https://xkcd.com/color/rgb/
                             ax=ax1)
            else:
                sns.lineplot(x=df.Position,
                            y=df.Predicted_curvature,
                            label="Curvature",
                            linewidth=1.9,
                            #color="C0",
                            color="xkcd:cerulean", #see 
                            # https://matplotlib.org/tutorials/colors/colors.html &
                            # https://xkcd.com/color/rgb/
                            ax=ax1)
                sns.lineplot(x=df.Position, 
                             y=df.Bendability,
                             label="Bendability",
                             linewidth=1.9,
                             #color="C1",
                             color="xkcd:orange", #see 
                            # https://matplotlib.org/tutorials/colors/colors.html & 
                            # https://xkcd.com/color/rgb/
                             ax=ax1)
            plt.legend()
        else:
            # Note that here I use the columns for the other reported metric
            # by referencing the end of the columns list because the first few
            # columns may be different depending if sequence was included.
            ax2 = ax1.twinx() # complexity and GC content need their own scale 
            # so plot on secondary y-axis; use of secondary y-axis with seaborn based on 
            # https://towardsdatascience.com/a-step-by-step-guide-for-creating-advanced-python-data-visualizations-with-seaborn-matplotlib-1579d6a1a7d0
            if smooth_plot_curves:
                # use of Akima1DInterpolator based on UPDATE#2 at 
                # https://stackoverflow.com/a/37223558/8508004 and seems to match
                # fairly well to Excel's 'Scatter with Smooth Lines' option
                import numpy as np
                from scipy.interpolate import Akima1DInterpolator
                akima1 = Akima1DInterpolator(df.Position, df.Predicted_curvature)
                akima2 = Akima1DInterpolator(df.Position, df[df.columns[-1]])
                x_4smooth = np.linspace(
                    min(df.Position), len(merged_sequence), 1000)
                sns.lineplot(x=x_4smooth,
                            y=akima1(x_4smooth),
                            label="Curvature",
                            linewidth=1.9,
                            #color="C0",
                            color="xkcd:cerulean", #see 
                            # https://matplotlib.org/tutorials/colors/colors.html &
                            # https://xkcd.com/color/rgb/
                            ax=ax1)
                sns.lineplot(x=x_4smooth, 
                             y=akima2(x_4smooth),
                             label=df.columns[-1],
                             linewidth=1.9,
                             color="xkcd:orange", #see 
                            # https://matplotlib.org/tutorials/colors/colors.html &
                            # https://xkcd.com/color/rgb/
                             ax=ax2)
            else:
                sns.lineplot(x=df.Position,
                            y=df.Predicted_curvature,
                            label="Curvature",
                            color="xkcd:cerulean",
                            ax=ax1)
                sns.lineplot(x=df.Position, 
                             y=df.columns[-1],
                             label=df.columns[-1],
                             color="xkcd:orange",
                             ax=ax2)
            # getting both legends as one based on http://matplotlib.1069221.n5.nabble.com/legend-and-twiny-td14497.html
            #handles1, labels1 = ax1.get_legend_handles_labels()
            #handles2, labels2 = ax2.get_legend_handles_labels()
            #ax2.legend(handles1+handles2, labels1+labels2)

            # getting both legends as one based on https://stackoverflow.com/a/14344146/8508004
            #h1, l1 = ax1.get_legend_handles_labels()
            #h2, l2 = ax2.get_legend_handles_labels()
            #ax1.legend(h1+h2, l1+l2, loc=2)

            # based on https://stackoverflow.com/questions/5484922/secondary-axis-with-twinx-how-to-add-to-legend
            #ax1.legend(loc=0)
            #ax2.legend(loc=1) 
            
            #legends for both lines made automatically, based on 
            # https://stackoverflow.com/a/47370214/8508004
            fig.legend(
                loc="upper right", bbox_to_anchor=(1,1), 
                bbox_transform=ax.transAxes)
            # While the gathered legends works, it seems when using seaborn
            # the legends from each of the two axes show up on the plot over each other 
            # still. Can use the code described at https://stackoverflow.com/a/54782540/8508004
            # along with the addition in the comment by Mihai Cherlau to hide 
            # those. (also see https://stackoverflow.com/a/38136752/8508004 )
            # Note: before deleting can use `dir(ax1.legend_)` to see all 
            # attributes and methods associated with the seaborn `legend_`.
            ax1.legend_.remove()
            ax2.legend_.remove()
        plt.ylabel("");

        right_side_buffer = 0
        plt.xlim(0, len(merged_sequence)+right_side_buffer)

        # add slashes, |, and parenthese back of special description line 
        # pattern recognized from the sanitized pattern. In other words, make 
        # the title text look like the actual input description line that 
        # includes date stamp in user special case example `R5|11/23/19|16(4)`,
        # which gets sanitized to `R5_11-23-19_16+4+`
        if show_date_with_slashes_in_plot_title and (
            sample_name.count("_") == 2) and (
            sample_name.split("_")[1].split("_")[0].count("-") == 2) and (
            sample_name.split("_",2)[2].count("+") == 2) and (
            sample_name.endswith("+")):
            title_text = sample_name.split(
                "_")[0]+"|" + sample_name.split("_")[1].replace(
                "-","/") +"|" + sample_name.split("_",2)[2]
            #restore parantheses near end too. Can use fact `+` at end to
            # to specify closing parantheses. So position of the `+` used to 
            # enable using `+` as a symbol for two different characters.
            if title_text.endswith("+") and "+" in title_text[:-1]:
                title_text = title_text[:-1] + ")"
                title_text = title_text.replace("+","(")
        else:
            title_text = sample_name


        #plt.title(f"{sample_name}",color="#333333",fontsize = 16)
        #plt.title(f"{sample_name}",color="#333333",fontsize = 'large', fontweight='bold')
        plt.title(f"{title_text}",color="#333333",fontsize = 15, 
            fontweight='bold')

        plt.show()
        sns.set() #resets to default, in case want to plot other things below; 
        # inherited from code for twin y-axes (secondary y-axis) with seaborn 

        # And starting another spinner back up now allows the previous plot
        # to remain as more data processed and sets up for next plot to display
        spinner = Halo(
            text=f'Processing file #{indxf+1} entry #{indxe+1} ...', 
            spinner='dots',color = 'magenta')
        spinner.start()


    


        # CLEAN UP:
        #----------------------------------------------------------------------#
        # Clean up by:
        # 1. deleting the temporarily made fasta file
        # 2. deleting the output files made, should be a png and text file
        # 3. reset log_file_text if lightweight
        #!rm {fasta_file_name_for_merged}
        #!rm {output_file_suffix}
        os.remove(fasta_file_name_for_merged)
        os.remove(output_file_suffix)
        if (not include_gnuplots) or (lightweight_archive):
            #!rm {output_file_suffix}.png
            os.remove(f"{output_file_suffix}.png")

        if lightweight_archive:
            log_file_text = lw_note_text + "\n\n    --LOGS OF INDIVIDUAL \
            SEQUENCE PROCESSING SCRUBBED TO KEEP LIGHTWEIGHT--\n\n\n"

        

        # SET-UP RESULTS FOR STORAGE LATER:
        #----------------------------------------------------------------------#
        #store df
        #store plot
        prefix_4_df_pkl_n_plot_file_saves = f"{sample_set_name}_{sample_name}"
        df_pkl_nom = f"{prefix_4_df_pkl_n_plot_file_saves}.pkl"
        df.to_pickle(df_pkl_nom)
        files_produced.append(df_pkl_nom)
        df_tsv_nom = f'{prefix_4_df_pkl_n_plot_file_saves}.tsv'
        df.to_csv(df_tsv_nom, sep='\t',index = False)
        files_produced.append(df_tsv_nom)

        # if lightweight_with_images used then want to make the images, but 
        # since lightweight still don't include them in the archive
        if make_images:
            plot_png_nom = f"{prefix_4_df_pkl_n_plot_file_saves}.png"
            fig.savefig(plot_png_nom)
            plot_svg_nom = f"{prefix_4_df_pkl_n_plot_file_saves}.svg"
            fig.savefig(plot_svg_nom) #vector graphics best for directly scaling
        if not lightweight_archive:
            files_produced.append(plot_png_nom)
            files_produced.append(plot_svg_nom)

        
        dataframes_produced_keyed_by_sample[sample_name] = df
        if lightweight_archive:
            plots_produced_keyed_by_sample[sample_name] = "saving serialized \
            Python plot object scrubbed due to `lightweight` setting"
        else:
            plots_produced_keyed_by_sample[sample_name] = fig

        if include_gnuplots and (not lightweight_archive):
            files_produced.append(f"{output_file_suffix}.png")


    spinner.stop()


    # PREPARE SAMPLE SET SUMMARY:
    #--------------------------------------------------------------------------#
    # Show this in the stderr output at end of run and also send to Log file.
    # As part of it, it includes some quality control checks & looks for 
    # duplicate sequences. Because mode calculation can produce errors if there 
    # is not actual value that is most frequent some special handling of that 
    # has been built in based on https://stackoverflow.com/a/25109040/8508004. 
    for_out2both = "\n\n#-----------******SUMMARY*******-----------#"
    log_file_text = out2_stderr_n_log(for_out2both,log_file_text)
    for_out2both = "\nSAMPLE SET:"
    log_file_text = out2_stderr_n_log(for_out2both,log_file_text)
    for_out2both = f"\n{sample_set_name}"
    log_file_text = out2_stderr_n_log(for_out2both,log_file_text)
    #----SUMMARIZE NUMBER & UNIQUENESS-----------------------------------------#
    for_out2both = "\n#-----*NUMBER & UNIQUENESS*-----#:"
    log_file_text = out2_stderr_n_log(for_out2both,log_file_text)
    processed_seqs=[x for x in sequences_processed_keyed_on_name.values()]
    len_processed_seqs=[len(
        x) for x in sequences_processed_keyed_on_name.values()]
    # Some quality control checks.
    # in anticpation of producing a summary, check for duplicates
    number_unique_sequences_processed = len(set(processed_seqs))
    if len(sequences_processed_keyed_on_name) != len(set(processed_seqs)):
        sequences_with_duplicates = True
        merged_sequences_with_duplicates = True
        # warn in std.err if noted, even though data numbers of uniques may 
        # be in summary, too
        for_out2both = "\n***WARNING: Duplicate sequences detected.***"
        log_file_text = out2_stderr_n_log(for_out2both,log_file_text)
        for_out2both = ("\n{} cassette sequences processed; {} "
            "unique.".format(len(sequences_processed_keyed_on_name),
            number_unique_sequences_processed))
        log_file_text = out2_stderr_n_log(for_out2both,log_file_text)
    else:
        for_out2both = (
            "\n{} unique cassette sequences processed.".format(
            number_unique_sequences_processed))
        log_file_text = out2_stderr_n_log(for_out2both,log_file_text)
    processed_merged_seqs=[str(
        x) for x in merged_sequences_processed_keyd_on_name.values()]
    len_processed_merged_seqs=[len(
        str(x)) for x in merged_sequences_processed_keyd_on_name.values()]
    assert len(set(processed_seqs)) == len(set(processed_merged_seqs)), ("The "
        "number of sequences of the cassette elements should be the same as "
        "the number of merged sequences.")
    #----SUMMARIZE LENGTH------------------------------------------------------#
    for_out2both = "\n#-----------*LENGTH*------------#:"
    log_file_text = out2_stderr_n_log(for_out2both,log_file_text)
    for_out2both = ("\nThe max cassette size:\n{} nts\nMin cassette "
        "size:\n{} nts".format(max(len_processed_seqs),min(len_processed_seqs)))
    log_file_text = out2_stderr_n_log(for_out2both,log_file_text)
    for_out2both = ("\nThe mean cassette size:\n{:.0f} nts\nMedian cassette "
        "size:\n{:.0f} nts".format(
        mean(len_processed_seqs),median(len_processed_seqs)))
    log_file_text = out2_stderr_n_log(for_out2both,log_file_text)
    try:
        lenmode= mode(len_processed_seqs)
        for_out2both=("\nThe mode cassette size:\n{:.0f} nts".format(lenmode))
    except StatisticsError:
        for_out2both=("\nThere is no mode for cassette size.")
        if len(set(len_processed_seqs)) == len(len_processed_seqs):
            for_out2both=(" All values are unique.")
    log_file_text = out2_stderr_n_log(for_out2both,log_file_text)
    
    for_out2both = ("\nThe max size for merged with defined flanking "
        "sequences:\n{} nts\nMin size for merged with defined flanking "
        "sequences:\n{} nts".format(max(len_processed_merged_seqs),
        min(len_processed_merged_seqs)))
    log_file_text = out2_stderr_n_log(for_out2both,log_file_text)
    for_out2both = ("\nThe mean size for merged with defined flanking "
        "sequences:\n{:.0f} nts\nMedian size for merged with defined flanking "
        "sequences:\n{:.0f} nts".format(
        mean(len_processed_merged_seqs), median(len_processed_merged_seqs)))
    log_file_text = out2_stderr_n_log(for_out2both,log_file_text)
    try:
        lenmmode= mode(len_processed_merged_seqs)
        for_out2both=("\nThe mode size or merged with defined flanking "
            "sequences:\n{:.0f} nts".format(lenmmode))
    except StatisticsError:
        for_out2both=("\nThere is no mode for size of merged.")
        if len(set(len_processed_merged_seqs)) == len(
            len_processed_merged_seqs):
            for_out2both=(" All values are unique.")
    log_file_text = out2_stderr_n_log(for_out2both,log_file_text)

    #----SUMMARIZE %G+C:-------------------------------------------------------#
    cassette_nt_counts = {}
    to_add_seq_to_df = []
    for sn,seq_str in sequences_processed_keyed_on_name.items():
        cassette_nt_counts[sn] = collections.Counter(seq_str.upper())
        assert len(cassette_nt_counts[sn]) <= 4, ("There should only be "
            f"G,A,T,or C in count of cassette nts of {sn}.")
        expect_nts = ["G","A","T","C"]
        if len(cassette_nt_counts[sn]) < 4:
            # add zeroes for any missing nucleotides so later uses of those 
            # columns don't cause errors
            for nt in expect_nts:
                if nt not in cassette_nt_counts[sn]:
                    cassette_nt_counts[sn][nt] = 0
        to_add_seq_to_df.append(seq_str.upper()) 
    cassette_nt_count_df = pd.DataFrame.from_dict(
        cassette_nt_counts, orient='index').fillna(0)
    # It seems for columns that include a zero fill on nt count, the dtype gets 
    # converted to a float, this will lock back into integer
    cassette_nt_count_df.G = cassette_nt_count_df.G.astype(dtype='int64')
    cassette_nt_count_df.A = cassette_nt_count_df.A.astype(dtype='int64')
    cassette_nt_count_df.C = cassette_nt_count_df.C.astype(dtype='int64')
    cassette_nt_count_df["T"] = cassette_nt_count_df["T"].astype(dtype='int64')#
    # for nucleotide `T` need to use bracket notation and not attribute notation
    # a.k.a. dot notation because `DataFrame.T` is something in Pandas for 
    # transposing, see https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.T.html
    
    # include the sequence for cassettes as they are probably short and move to
    # after the name
    cassette_nt_count_df["seq"] = to_add_seq_to_df 
    # move `seq` column to first in dataframe based on 
    # https://stackoverflow.com/a/51009742/8508004
    cols = cassette_nt_count_df.columns.tolist()
    n = int(cols.index('seq'))
    cols = [cols[n]] + cols[:n] + cols[n+1:]
    cassette_nt_count_df = cassette_nt_count_df[cols]
    cassette_nt_count_df["Total_nts"] = cassette_nt_count_df.sum(1)
    cassette_nt_count_df['%G+C'] = (
        cassette_nt_count_df[['C','G','Total_nts']].apply(
        percent_GCcalc, axis=1))
    cassette_nt_count_df["rank_on_GC"]= cassette_nt_count_df['%G+C'].rank()#from 
    # `GSD Asessing_ambiguous_nts_and_calculating_percentGC_for_332_saccharomytina_genomes.ipynb`
    cassetteGClist = cassette_nt_count_df['%G+C'].tolist()
    merged_nt_counts = {}
    for sn,seq_str in merged_sequences_processed_keyd_on_name.items():
        merged_nt_counts[sn] = collections.Counter(seq_str.upper())
        assert len(merged_nt_counts[sn]) <= 4, ("There should only be "
            f"G,A,T,or C in count of nts of merged {sn}.")
        if len(merged_nt_counts[sn]) < 4:
            # add zeroes for any missing nucleotides so later uses of those 
            # columns don't cause errors
            for nt in expect_nts:
                if nt not in merged_nt_counts[sn]:
                    merged_nt_counts[sn][nt] = 0
    merged_nt_count_df = pd.DataFrame.from_dict(
        merged_nt_counts, orient='index').fillna(0)
    # It seems for columns that include a zero fill on nt count, the dtype gets 
    # converted to a float, this will lock back into integer. Less likely to
    # happen for merged form but might as well avoid the chance.
    merged_nt_count_df.G = merged_nt_count_df.G.astype(dtype='int64')
    merged_nt_count_df.A = merged_nt_count_df.A.astype(dtype='int64')
    merged_nt_count_df.C = merged_nt_count_df.C.astype(dtype='int64')
    merged_nt_count_df["T"] = merged_nt_count_df["T"].astype(dtype='int64')#
    # for nucleotide `T` need to use bracket notation and not attribute notation
    # a.k.a. dot notation because `DataFrame.T` is something in Pandas for 
    # transposing, see https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.T.html
    merged_nt_count_df["Total_nts"] = merged_nt_count_df.sum(1)
    merged_nt_count_df['%G+C'] = (
        merged_nt_count_df[['C','G','Total_nts']].apply(
        percent_GCcalc, axis=1))
    merged_nt_count_df["rank_on_GC"]= merged_nt_count_df['%G+C'].rank()#from 
    # `GSD Asessing_ambiguous_nts_and_calculating_percentGC_for_332_saccharomytina_genomes.ipynb`
    # store the GC dataframes as files
    prefix_4_GCdf_file_saves = f"{sample_set_name}_cassettesGC"
    GCdf_pkl_nom = f"{prefix_4_GCdf_file_saves}.pkl"
    cassette_nt_count_df.to_pickle(GCdf_pkl_nom)
    files_produced.append(GCdf_pkl_nom)
    GCdf_tsv_nom = f'{prefix_4_GCdf_file_saves }.tsv'
    cassette_nt_count_df.to_csv(GCdf_tsv_nom, sep='\t')
    files_produced.append(GCdf_tsv_nom)
    prefix_4_mGCdf_file_saves = f"{sample_set_name}_mergedGC"
    mGCdf_pkl_nom = f"{prefix_4_mGCdf_file_saves}.pkl"
    merged_nt_count_df.to_pickle(mGCdf_pkl_nom)
    files_produced.append(mGCdf_pkl_nom)
    mGCdf_tsv_nom = f'{prefix_4_mGCdf_file_saves }.tsv'
    merged_nt_count_df.to_csv(mGCdf_tsv_nom, sep='\t')
    files_produced.append(mGCdf_tsv_nom)
    for_out2both = "\n#------------*%G+C*-------------#:"
    mergedGClist = merged_nt_count_df['%G+C'].tolist()
    log_file_text = out2_stderr_n_log(for_out2both,log_file_text)
    for_out2both = ("\nThe max cassette %G+C:\n{:.2%}\nMin cassette "
        "%G+C:\n{:.2%} ".format(max(cassetteGClist),min(cassetteGClist)))
    log_file_text = out2_stderr_n_log(for_out2both,log_file_text)
    for_out2both = ("\nThe mean cassette %G+C:\n{:.2%}\nMedian cassette "
        "%G+C:\n{:.2%}".format(mean(cassetteGClist),median(cassetteGClist)))
    log_file_text = out2_stderr_n_log(for_out2both,log_file_text)
    # round cassetteGClist entries for better chance at getting a mode? Note b/c
    # G+C is a percent need to shift rounding over two more decimal places
    rnd_cassetteGClist = [round(x, 3) for x in cassetteGClist]
    try:
        gcmode= mode(rnd_cassetteGClist)
        for_out2both=("\nThe mode rounded %G+C for cassettes"
            ":\n{:.1%}".format(gcmode))
    except StatisticsError:
        for_out2both=("\nThere is no mode for rounded %G+C of cassettes.")
        if len(set(rnd_cassetteGClist)) == len(rnd_cassetteGClist):
            for_out2both+=(" All values are unique.")
    log_file_text = out2_stderr_n_log(for_out2both,log_file_text)

    for_out2both = ("\nThe max %G+C for merged with defined flanking "
        "sequences:\n{:.2%}\nMin %G+C for merged with defined flanking "
        "sequences:\n{:.2%}".format(max(mergedGClist),min(mergedGClist)))
    log_file_text = out2_stderr_n_log(for_out2both,log_file_text)
    for_out2both = ("\nThe mean %G+C for merged with defined flanking "
        "sequences:\n{:.2%} nts\nMedian %G+C for merged with defined flanking "
        "sequences:\n{:.2%} nts".format(
        mean(mergedGClist), median(mergedGClist)))
    log_file_text = out2_stderr_n_log(for_out2both,log_file_text)
    # round mergedGClist entries for better chance at getting a mode? Note b/c
    # G+C is a percent need to shift rounding over two more decimal places
    rnd_mergedGClist = [round(x, 3) for x in mergedGClist]
    try:
        mgcmode= mode(rnd_mergedGClist)
        for_out2both=("\nThe mode rounded %G+C for merged with defined "
            "flanking sequences:\n{:.1%}".format(mgcmode))
    except StatisticsError:
        for_out2both=("\nThere is no mode for rounded %G+C of merged.")
        if len(set(rnd_mergedGClist)) == len(rnd_mergedGClist):
            for_out2both+=(" All values are unique.")
    log_file_text = out2_stderr_n_log(for_out2both,log_file_text)

    for_out2both = "\n#--------*****END OF SUMMARY******---------#"
    log_file_text = out2_stderr_n_log(for_out2both,log_file_text)
    



    # PREPARE DATA COLLECTION FROM THE SAMPLE SET FOR INCLUSION IN ARCHIVE:
    #--------------------------------------------------------------------------#

    seqs_dfs_and_plots_per_sample_set[sample_set_name] = [
        sequences_processed_keyed_on_name, 
        merged_sequences_processed_keyd_on_name,
        cassette_nt_count_df,
        merged_nt_count_df, 
        dataframes_produced_keyed_by_sample,
        plots_produced_keyed_by_sample
        ]

    #plus delete the `.fai` files made when processing current sequence file
    fai_fns = glob.glob('*.fai')
    [os.remove(x) for x in fai_fns]





# PACKAGE UP LOG AND ARCHIVE RESULTS FOR ALL:
#------------------------------------------------------------------------------#

now = datetime.datetime.now()

#store `sequences_processed_keyed_on_name` somehow
#store `merged_sequences_processed_keyd_on_name` somehow

# collect all the unsanitized input files if there are any
files_pattern_at_start_of_name = "unsanitized_input_"
unsanitized_input_fns= []
for file in os.listdir('.'):
    if fnmatch.fnmatch(file, files_pattern_at_start_of_name+'*'):
        unsanitized_input_fns.append(file)
if unsanitized_input_fns:
    files_produced.extend(unsanitized_input_fns)



#Serialize the data in `seqs_dfs_and_plots_per_sample_set`
# use of `with` based on https://stackoverflow.com/a/20101064/8508004
import pickle 
serial_fn = "seqs_dfs_and_plots_for_each_set.pkl"
with open(serial_fn, "wb") as f:
    pickle.dump(seqs_dfs_and_plots_per_sample_set, f)
files_produced.append(serial_fn)
if lightweight_archive:
    for_out2both = ("\n\nWhat was procesed and results have\nbeen "
        f"serialized for later use as `{serial_fn}`.")
else:
    for_out2both = ("\n\nIn addition to individal tabular data and plot "
        "images, records of\nwhat was procesed and results have\nbeen "
        f"serialized for later use as `{serial_fn}`.")
log_file_text = out2_stderr_n_log(for_out2both,log_file_text)


# Save the current notebook and store it as part of the archive next. <==DOESN'T
# WORK IN JUPYTERLAB, see https://github.com/jupyterlab/jupyterlab/issues/7627
#save_notebook("index.ipynb")
#files_produced.append("index.ipynb")
# Since cannot save the notebook, in JupyterLab, send the plots to a notebook
# using Jupytext and then execute & save that notebook with the archive?
if not lightweight_archive:
    plots4review_fn = make_and_run_review_nb(now,review_nb_stub, serial_fn)
    files_produced.append(plots4review_fn)
    # next line removes the `.py` intermediate made in process to make 
    # 'review' nb
    #!rm {plots4review_fn[:-6]+".py"}
    os.remove(plots4review_fn[:-6]+".py")
    for_out2both = ("\nA Jupyter notebook listing the resulting plots for "
        f"convenient reviewing\nhas been saved as `{plots4review_fn}`.")
    log_file_text = out2_stderr_n_log(for_out2both,log_file_text)


log_file_text += (f"\n\n#### THIS CONCLUDES THE RUN OF bendit analysis "
    f"{now.strftime('%b%d%Y%H%M')}#######")
log_file_name = f"LOG_ba{now.strftime('%b%d%Y%H%M')}.txt"
write_string_to_file(log_file_text, log_file_name )
sys.stderr.write(f"\n\nLog file from the run has been saved as "
    f"`{log_file_name}`.")
files_produced.append(log_file_name)

# Make the archive from the run
archive_file_name = f"bendit_analysis{now.strftime('%b%d%Y%H%M')}.tar.gz"
#print(files_produced) # FOR DEBUGGING ONLY
#!tar czf {archive_file_name} {" ".join(files_produced)}
os.system(f'tar czf {archive_file_name} {" ".join(files_produced)}')
sys.stderr.write("\n*****************DONE***********************************\n"
    "{} generated. Download it.\n"
    "*****************DONE***********************************".format(
    archive_file_name))


# CLEAN UP SO DIRECTORY NOT SO CLUTTERED SO MORE OBVIOUS WHAT TO DOWNLOAD:
#------------------------------------------------------------------------------#  
# Delete the files that have been archived (except leave the demo sequence file 
# so that development example can be run again easily)
if cleaning_step:
    sys.stderr.write("\nCleaning up...be patient...\n")
    spinner = Halo(text='Cleaning up...', spinner='dots',color = 'magenta')
    spinner.start()  
    if demo_file_name in files_produced:
        files_produced.remove(demo_file_name)
    for ef in files_produced:
        #!rm {ef}
        os.remove(ef)
    # Delete the `temp` files from the making of the dataframes step; they were used 
    # for the checking of the two sets of data bendIt reports.
    #!rm results_parsing_temp*.tsv
    results_parsing_temp_tsvs = glob.glob('results_parsing_temp*.tsv')
    [os.remove(rtsv) for rtsv in results_parsing_temp_tsvs]

    spinner.stop()
    sys.stderr.write("\nCleaning complete.")
#######------------------END OF MAIN SECTION------------------------------######
################################################################################
