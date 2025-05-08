# ReadMe Manual

## Table of Contents

- 1.0 General Information
   - 1.1 License Information
   - 1.2 Citing this Work
   - 1.3 Operating Google Colab
   - 1.4 Syntax Overview
   - 1.5 Graphical User Interface
- 2.0 Main Program
   - Step 1: Curation of Protein List using PDB
   - Step 2: Converting Accession Codes to Ensembl IDs
   - Step 3: Generating the Unique Gene List
   - Step 4: Retrieval of CPDB Databank and Other Variables
   - Step 5: Overrepresentation Analysis
   - Step 6: Hierarchical Clustering and Cluster Analysis Output
- 3.0 Secondary Program
   - Step 7: Downloading Specific Sets
- 4.0 Troubleshooting
   - 4.1 Changing File Names
   - 4.2 API and Web Page Issues


## 1. General Information

This manual provides detailed instructions on the use of the program titled “Cofactor Mapping &
Analysis Program (CoMAP).ipynb”.

The purpose of this code is to search the Protein Data Bank (PDB) for molecule entries that meet
the user specified criteria, generate a list of their unique Ensembl gene IDs, perform an
overrepresentation analysis on the biological pathways from ConsensusPathDB (CPDB) using
these gene IDs, and finally, cluster these significant pathways based on similarities between gene
sets. This program provides a novel way to study cofactor-dependent molecules.

This code has been created and run in the Google Colab environment using the Python 3 language.

Note: ChatGPT with GPT-4o was used as a tool to aid in the development of the Python script.

### 1.1 License Information

The Cofactor Mapping & Analysis Program (CoMAP) is freely available for academic, non-commercial use. Commercial users should contact Dr. Travis Craddock at
travis.craddock@uwaterloo.ca. Data used by CoMAP is available under the license terms of each
contributing databank (the Protein Data Bank, DAVID, and ConsensusPathDB).

### 1.2 Citing this Work

If using CoMAP in your research, please contact Dr. Travis Craddock at
travis.craddock@uwaterloo.ca for citation details. A publication is currently in preparation.

### 1.3 Operating Google Colab

To run a cell in Google Colab, click the start button located in the upper left-hand corner of the cell.
After clicking start, this button becomes a stop button.

A successfully run cell will display a green check mark next to the start button.

An unsuccessfully run cell will display a red exclamation mark next to the start button. This will be
accompanied by an error message in the terminal at the bottom of the code.

For this code specifically, the code blocks are sequential. This means that the first block of code,
the Main Program, must be successfully run before the Secondary Program can be.

### 1.4 Syntax Overview

This section provides a brief overview of some important syntax for this code.

Lines of text beginning with a hashtag (#) are code comments. Code comments are not run when
the code is executed and do not alter the execution of the code. These comments are intended to
provide information on what specific lines or sections of code are doing.

Any line of code following the format def function_name(input 1, input2, etc.): is used to define a
function. These functions are what complete the tasks necessary to run the code.


Any line of code following the format variable_name = “File name.extension” is used to define an
output file. The file name portion of these lines are the only pieces of code that should be changed,
though it is not recommended.

Please read Section 4.1 before altering file names.

### 1.5 Graphical User Interface

The Main Program and Secondary Program both have a graphical user interface (GUI). The GUI will
be output in the terminal at the bottom of each block of code upon successful execution.

For the Main Program GUI, the first dropdown selects the source organism, the organism which the
molecule originally came from, with the options being: Homo sapiens (human), Mus musculus
(mouse), and Saccharomyces cerevisiae (yeast). Next, you can select the polymer entity type. A
default entity is already present and cannot be deleted. The dropdown menu defaults to “protein”,
but you may change it to any of the supported PDB polymer types: protein, DNA, RNA, nucleic acid
hybrid (NA-hybrid), or other. The green “Add Polymer Entity” button and red “Remove Polymer
Entity" button allow you to submit between 1 and 5 polymer types. The grey “AND/OR” button
determines how multiple polymer types are interpreted with “AND” querying structures that contain
all selected polymer types and “OR” querying structures that contain any of the selected polymer
types. The next fields are for cofactor selection. One cofactor entry is automatically created and
must be submitted. The cofactor name can be chosen by either clicking “Dropdown” and using one
of the pre-selected options or clicking "Custom” and entering the desired cofactor name into a
textbox. Custom entry text must match with cofactor names available in the PDB. More cofactors
can be added using “Add Cofactor” and removed using “Remove Cofactor”. The “AND/OR Button”
dictates how two or more cofactors will be combined in the search. An “AND” selection will search
for molecules that contain all the specified cofactors while an “OR” selection will look for
molecules that contain at least one of the cofactors. The threshold q-value field is a textbox that
takes in a decimal number. A decimal point or scientific notation (using E or e) can be used to enter
a number. A zero or negative value will change to 0.05 after submission. The last four true or false
dropdowns are used to decide whether files from various steps in the code will be downloaded.

Some important notes about this GUI are as follows. The first time you run the Main Program in a
session several packages and modules will need to be downloaded. The program will do this
automatically. The terminal will update you on the progress of the installations and downloads, but
this process can take up to a minute. During this process, however, nothing is downloaded or
installed onto your physical computer. The only information that will be downloaded to your
computer are the selected files during program execution. Additionally, the “Submit” button will
become disabled after it is pressed. This means values cannot be changed and resubmitted.

For the Secondary Program GUI, the first textbox takes in the cutoff value as a decimal number. This
can again be entered by typing out the decimal number or using scientific notation. The second
textbox takes in an integer number for the minimum number of pathways needed in a group to be a
viable cluster. A negative or zero number submission will update to 0.1 for the cutoff value and 1 for
the pathways value upon submission. The first true or false dropdown is for whether or not to show
the dendrogram. The dendrogram is clearer and more legible when there are fewer overrepresented
pathways. The show gene cluster membership field dictates whether or not to print the Ensembl ID


of the genes and the clusters they are a part of. Additional input fields that appear after selecting
true on this option allow you to specify which genes will be printed to the terminal, depending on
how many clusters the gene appears in. The last dropdown is for downloading the set file.

There are several important things to know about the execution of this GUI. The Secondary Program
can only be successfully run after the Main Program has finished running and produced a heat map.
If too few pathways are overrepresented or a heat map is not produced, this code won’t execute.
The Secondary Program must also be run in the same session as the finished execution of the Main
Program. If the program times out and disconnects, the Secondary Program will not have the
information it needs to execute. Finally, this code can be run multiple times. The submit button for
this GUI does not automatically disable meaning that multiple values can be submitted, allowing
you to download all the information you require.

## 2. Main Program:

### Step 1: Curation of Protein List using PDB

The search query used to filter through the Protein Data Bank is created. It consists of the source
organism, polymer entity type(s), cofactor(s), and the and/or operators. All entries meeting the
query criteria are returned and saved in a single file. This file contains information on the entry ID,
PDB ID, gene name, macromolecule name, and accession code(s) of each entry. This file is referred
to as the PDB entries file in the GUI.

### Step 2: Converting Accession Codes to Ensembl IDs

The accession codes from the PDB entries file are submitted to the DAVID Gene ID Conversion Tool.
The Ensembl gene IDs for each accession code are appended to and saved as a single txt file.

### Step 3: Generating the Unique Gene List

The PDB entries file will be updated to include a new column for the Ensembl gene IDs. The
Ensembl gene IDs will match up with their corresponding accession code. A list of the unique
Ensembl gene IDs retrieved will be generated. The gene names will be arranged in alphabetical
order from uppercase A to Z followed by lowercase a to z. The updated PDB entries file and unique
Ensembl gene list can be downloaded.

### Step 4: Retrieval of CPDB Databank and Other Variables

The code will automatically navigate to the ConsensusPathDB organism databank for your selected
source organism. All the biological pathways contained in CPDB for that organism will be
downloaded with their genes identifiable by Ensembl ID. This file is initially a tab file, which gets
converted to a csv file.

The program also retrieves two numeric values from CPDB. The first is the total background size (N).
This variable represents the total number of Ensembl IDs that are in at least one CPDB pathway.
The second is the number of mapped entities (M). This number signifies the number of genes from
the unique gene list that are present in at least one CPDB pathway.


Please read Section 4.2 API and Web Page Issues if you are encountering problems with CPDB.

### Step 5: Overrepresentation Analysis

An overrepresentation analysis is run using the unique Ensembl ID list. The match percentage (the
percentage of genes in a pathway that are found in the unique gene list), the p-value (the measure
of how likely a result is due to chance rather than being a meaningful biological result), and the q-value (the measure of if a result is still significant after adjusting for multiple tests) are calculated.
The pathways are then filtered according to their q-value, with any pathway having a q-value greater
than the threshold value being deemed insignificant and being discarded. The remaining pathways
are added to a new file, the filtered p- and q-value file. This file includes columns on the pathway
names, pathway sources, pathway IDs, Ensembl IDs from the unique gene list that are a part of the
pathways, the match percentages, and the p- and q-values. This file can be downloaded.

Note: The q-values generated by the code may differ slightly (by an order of magnitude or two) from
those generated by ConsenusPathDB. The p-values, however, should not differ significantly, if at all.

### Step 6: Hierarchical Clustering and Cluster Analysis Output

The overrepresented pathways are clustered together based on the cutoff value and minimum
number of pathways value. The cutoff value is how dissimilar two pathways are allowed to be from
each other. The lower the value, the more similar the pathways’ gene sets need to be for them to be
clustered together. The minimum number of pathways value specifies how many pathways need to
be in a group for it to be considered a cluster. A unique combination of these two variables will be
called a set. A total of 152 sets are used in the heat map. The minimum number of pathways value
ranges from 1 to 8 in whole number increments and the cutoff value ranges from 0.05 to 0.95 in 0.
increments. A final cluster analysis file can de downloaded. This file will contain information on all
152 sets. For every set the file will have information on its number of clusters, percentage of
pathways expressed, percentage of genes expressed, average compactness, and word frequency.
The average compactness of a set is a measure of how closely related the pathways in each cluster
are. The word frequency is the single most frequent word(s), excluding common filler words or other
unhelpful descriptive terms, used in the pathway names of each cluster.

## 3. Secondary Program

### Step 7: Downloading Specific Sets

The set downloaded depends on the cutoff value and minimum number of pathways value you
input in the GUI. This downloadable file contains columns on the cluster IDs, the pathways in the
clusters, the compactness score of the clusters, and the cluster labels.

Note: This program needs to be run directly after the Main Program and in the same session so that
the functions can properly take in the necessary variables. New values can be submitted as many
times as you desire after the Main Program has finished executing.


## 4. Troubleshooting

### 4.1 Changing File Names

Changing file names may be necessary. It is recommended that you change file names after they
have been downloaded to your computer as opposed to altering the code. If, however, it is
necessary for you to modify how files are named in the code, below is some crucial information to
ensure that this is done properly.

The naming of a file in this program has three parts: 1 - the variable name, 2 - the file name, and 3 -
the extension. Do not alter variable names or extensions. Variable names are used in multiple
places in the code and changing one instance can have repercussions on the functioning of the
program. All file extensions are csv. The csv module is the one that has been used in the creation of
this program and changing this can also affect code function. Csv files should be openable in Excel
and can then be saved with different extensions. Therefore, only the file name section should be
altered.

Some file names contain variables within them. Variables appear within curly brackets ({ }). Do not
change the variables within the curly brackets. However, you can delete the variables in their
entirety (including the curly brackets) and change the file name in that manner.

If a file name has “/tmp/” at the beginning, do not modify that portion of the name. This indicates
that the file is located in a temporary download folder and must be accessed from there for proper
use and download.

Text that does not meet any of the above-mentioned criteria can be safely changed.

The five variable names that correspond to the five downloadable files—and the only variables that
should be changed—are: combined_report, output_csv_path, filtered_P_and_Q_file, cluster_file, and
specific_path_cut_file.

### 4.2 API and Web Page Issues

Something to be aware of regarding this program is that the external databanks (PDB, DAVID, and
CPDB) are accessed every time during code execution. Due to this, these databanks must be
functioning at the time the code is being run. If one of the websites is down or unable to be
accessed for another reason, the code will not execute.

You can check if one of the required webpages are down at these links:

- The Protein Data Bank: https://www.rcsb.org/
- DAVID Conversion Tool: https://davidbioinformatics.nih.gov/conversion.jsp
- Consensus Path DB – Human: http://cpdb.molgen.mpg.de/
- Consensus Path DB – Yeast: http://cpdb.molgen.mpg.de/YCPDB
- Consensus Path DB – Mouse: http://cpdb.molgen.mpg.de/MCPDB

The human databank at CPDB can sometimes be down for multiple days at a time. In that case, the
human organism will not work with the program, though the other organisms may still be available.


Unfortunately, in the case that one of the above websites is down for maintenance, has reached
server capacity, or is otherwise unable to be accessed, the only thing that can be done is to wait for
proper functioning to be restored then try again.


