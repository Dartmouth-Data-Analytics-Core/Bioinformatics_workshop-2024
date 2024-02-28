# Unix/Linux Shell basics

The Unix/Linux *'Shell'* describes a program that takes commands from an input (eg. your keyboard) and passes them to an operating system that will execute them. In contrast to a *Graphical User Interface (GUI)* the Shell is simultaneously a *command line interface (CLI)* and a programming language that allows you to perform tasks on your system.

<p align="center">
  <img src="../figures/terminal.png" height="75%" width="75%"/>
</p>

Interacting with a system through the Shell has many advantages over a GUI. The Shell allows you to quickly and easily navigate through directories on your computer, make, copy and search files in a systematic way, and construct pipelines that will execute complex tasks on big datasets.

Importantly, the Shell allows us to do each of these in the context of Bioinformatics, and Bioinformatics software.

## Why learn to use a Shell?  
Learning to use a Shell can be challenging, however it is a key skill in bioinformatics, as it is the primary way in which we interface with a lot of bioinformatics software and file types.

Some bioinformatics software provides GUIs that enable users to execute tasks with programs that you would otherwise execute using the Shell. While such software can be powerful in the right context, they can also make it very easy to perform tasks in bioinformatics incorrectly, and should therefore be treated with caution.


## The Bash shell

### The absolute basics

There are different types of Unix shells, however the most popular is Bash (the *Bourne Again Shell*), which is also the most common on Linux systems. Since the majority of participants will be using the Bash shell, and this is the default shell used on Dartmouth's high performance computing system (which we will be using), this lesson will introduce the Shell through using the Bash shell, however most, if not all, content should be transferable to other Unix shells.

> Use the Cheat Sheet in the GitHub repo to help you learn commands and available options.

Accessing the (bash) shell:  
- On a Mac or Linux system, the *Terminal* application provides access to the shell. There are also applications that you can download that provide customizations not present in the Terminal application, such as [iTerm2](https://iterm2.com/).
- On a Windows system, you can use an application such as [MobaXterm](https://mobaxterm.mobatek.net/) or [PuTTY](https://www.putty.org/).

<p align="center">
  <img src="../figures/shell.png" height="80%" width="80%"/>
</p>

When you open your terminal application you will be presented with the command prompt `$` when you are able to input commands. If the terminal is busy and cannot currently accept new commands, you will not be presented with the prompt.

When the prompt is shown, you can enter commands by typing them in after the prompt. Commands are typically composed of three components:  
- the name of the command itself  
- any flags or options you wish to run the command with (not always required)
- a file or directory to act on (sometimes implicit)

In the above example, we are asking the Shell to pass the `mkdir` command to the operating system (for making directories) with the `-p` option (which lets us make parent and sub directories at the same time) and the argument detailing what directory we want the command to make.

Manual pages for specific commands can be accessed using the `man` command.
```bash
man mkdir
```

The shell provides us with commands that allow us to list files in our current working directory, as well as change the current working directory to another location. For example:
```bash
# 'ls' command lists files in our current working directory
ls

# run ls with the '-a' option to include hidden files
ls -a

# The pwd command shows you your current working directory
pwd

# cd allows you to change your current working directory ('.' means current directory)
cd .

# '..' tells the shell to move your current directory up one directory
cd ..

# check you directory again
pwd

# now go back down the tree.  Replace the directory name with your own.
cd OwenW/
pwd
```

To go back down the directory structure, we specified a directory that was in our current working directory (cd). This is called a **relative path**, since it is relative to our current directory and will only work if our current directory is relative to the directory we are trying to reach.  

Relative paths are contrasted to **absolute paths** which always starts with a '/' and will start at the root (highest level) of the directory tree, and work from wherever you are in the directory substructure. For example:
```bash
ls /Users/OwenW/
```

By default, your terminal application will start your current directory as your *home directory* (more on that later). No matter where you are, you can always get back to your home directory using the tilde `~` with the `cd` command.
```bash
cd ~
```

Another useful command is `echo` which will evaluate and print characters provided to it.
```bash
echo "words words words"
```

We can use the redirect command (>) to redirect the output of commands like echo into a file. As an example, lets save the important note we made above to a text file.
```bash
echo "words words words" > mynotes.txt
```

## Log on to discovery cluster

Many of the higher level commands for working with NGS data will require a lot of memory and computing power, more than most laptops can handle efficiently.
The discovery cluster is a resource hosted by Dartmouth's Research Computing team. This cluster enables you to execute high level commands without using the memory and computing power on your local machine (more on this soon). Let's log onto the discovery cluster now. We will use a secure shell command `ssh` to log onto the discovery cluster.

```bash
# Establish the secure shell connection
ssh netID@discovery7.dartmouth.edu

# Enter your password at the prompt (when you type no characters will show up to preserve privacy)
netID@discovery7.dartmouth.edu's password:

# You're in!
(base) [netID@discovery7 ~]$

```
The commands that you just executed locally in your terminal window work the same way when you are logged into discovery. It is always useful to orient yourself when you're working on an HPC so that you know where the output of all of the commands you run will end up. Let's run our first command to get your location.

```bash
# Check your location on the cluster
pwd
```
You should see something like `/dartfs-hpc/rc/home/h/netID` displayed in response to your command. Initially when you log on you will always be directed to your home directory (the address or path listed above). Your home directory by default will have 50GB of storage space to begin with, if you are running something that requires more storage space it is possible to extend that limit temporarily with the `/dartfs-hpc/scratch/ drive`. This is where we have stored all of the files you will be working with today. Directories and files hosted on the `/dartfs-hpc/scratch/` drive will only be kept for 45 days, you will receive a notification from Research Computing before the data is deleted.


## Start a new directory

Let's start by making a folder, or directory, to store all of the work you do today. It's often a useful practice to give directories informative names, for example `230601_YOURINITIALS_rnawrkshp`, (ex. Owen's `230601_omw_rnawrkshp`, Shannon's `230601_sms_rnawrkshp`, Tim's `230601_sullivan_rnawrkshp`). At a glance I could tell that this directory contains data from June 01, 2023 pertaining to an RNAseq workshop. Notice I chose a title with no spaces. The 'space' is a *special character* in UNIX, and special characters need to be *escaped* with the `\` and so `230601_abc_workshp` would look like `230601\ abc\ workshp` with escape characters.

File names with spaces become unwieldy to type out so most programmers will replace spaces with `_`, `.`, or `-` in their filenames to keep everything neat.

```bash
# Navigate to scratch so you can make your own directory there
cd /dartfs-hpc/scratch/

# Make the directory.  Replace 'omw' with your own username.
mkdir 240313_sms_workshp/

# Change to the newly-created directory.
# replace omw with your initials
cd 240313_sms_workshp

# Check your location on the cluster
pwd

# List the contents of your directory
ls
```
As expected, the new directory that you created is empty there are no files. 

## Customize your environment

What is an alias
What is a .bash_profile

This can be a real time saver to save alias for frequently accessed work directories.

Add the following definitions to your `.bash_profile` using the `nano` text editor.
```
nano ~/.bash_profile
```
This opens your `.bash_profile` where you can copy the following definitions and add them to the bottom of the file, do not forget to edit the path of the **FOB** definition to include your own initials. 

```
# define the location for the workshop resources we provide
RESOURCE="/dartfs/rc/nosnapshots/G/GMBSR_refs/fobwrkshp/"

# define the location where you will be working during the workshop
FOB="/dartfs-hpc/scratch/240313_sms_fob"
```

Next use the `ctrl + X` keys to close the file, type `Y` to save your changes. Now let's reload your `.bash_profile` to customize your session and check that our definitions are working as expected.

```
source ~/.bash_profile

echo $FOB
echo $RESOURCE
```

### Viewing the contents of files

Lets copy a file from the resource directory we created for this workshop to the directory you just created. This file (`all_counts.txt`) provides raw read counts for an RNA-seq experiment, with genes in rows and samples in columns.

```bash
#it is prudent to run the 'pwd' command before running the 'cp' command so you know where your files will be copied to
pwd

# copy the file from the resource directory to the $FOB directory you just created

cp $RESOURCE/all_counts.txt $FOB

```

The shell provides us with commands to view the contents of files in defined ways. The `cat` command for example (which stands for for concatenate) will print the entire contents of a file to the terminal. This can be useful for smaller files, but as you will see with larger files can quickly fill the terminal with more lines of data than it can display.

```bash
cat all_counts.txt
```

When working with larger files, which we are usually doing in bioinformatics, you may not wish to print the whole file as it would overrun your terminal. Other commands exist that allow you to explore file contents with more control.
- `more` shows you as much of the file as can be shown in the size of the terminal screen you have open, and you can continue to "scroll" through the rest of the file by using the space bar  
- `less` is a similar command to `more`, and has advantages such as not persisting in the terminal, and being searchable
- `head` will print the first 10 lines by default, but this number can be controlled with the `-n` option
- `tail` will print the final 10 lines of a file, and can also be controlled with the `-n` option

We will use a large text file to show the utility of these commands, as well as other commands in the subsequent parts of this lesson.
```bash
# Show the first 10 lines of the all_counts.txt file
head -n 10 all_counts.txt

# Show the last 20 lines of the all_counts.txt file
tail -n 20 all_counts.txt

# Use the word count (wc) command with the lines option (-l) to show how many lines (rows) are in the dataset
wc -l all_counts.txt
```

### Renaming and removing files

Sometimes you will need to reorganize your directories or rename a file, which can be achieved with the `mv` command. Let's start by copying the all_counts.txt file from the rnaseq_workshp directory to your home directory.

```bash
# Copy the all_counts.txt file to your home directory
cp all_counts.txt ~/all_counts.txt
```
Now let's rename the copy of the all_counts.txt file that we just created.
```bash
# Rename the copied all_counts.txt file
mv ~/all_counts.txt ~/all_counts.copy.txt
```

You can also use the `mv` command to move a file to a new location. Let's move the all_counts.copy.txt from your home directory into your rnaseq_workshp directory.

```bash
# Move the all_counts.copy.txt into your rnaseq_workshp directory
# use pwd to check that you are in your own directory first
mv ~/all_counts.copy.txt ./

#check the contents of your rnaseq_workshp directory
ls
```

Copying the all_counts.copy.txt file was just an exercise to show you how the tools work, in practice you will want to keep your directories as neat as possible as you accumulate a lot of files. Let's remove the all_counts.copy.txt file with the `rm` command.

```bash
# For the sake of being careful, let's first list the details file to be removed
ls -l all_counts.copy.txt
# Remove the all_counts.copy.txt file
rm all_counts.copy.txt
```

You will notice that before the file was deleted you were asked if you were sure you wanted this file deleted. You want to be careful not to remove files that you did not create if you are working in shared directories. If you want to bypass this checkpoint, you can use the `-f` flag with `rm -f` to force the removal of a file, but be careful with this, as there is no *Trash* equivalent in the shell.

### Manipulating file contents

Some commands enable you to manipulate and subset files based on specific parameters. One useful example is the `cut` command, which allows you to 'cut' a file based on the options you select, such as the `-f` option, which corresponds to fields (columns). We could use `cut` to obtain read counts for only the first 5 samples in `all_counts.txt`.
```bash
# Look at only the counts from the first five columns
cut -f 1,2,3,4,5 all_counts.txt
```

To prevent all rows being printed to our console, we could combine the `cut` command with the `head` command using a *'pipe'*, specified by a '|'. Pipes send the output an initial command to a subsequent command, all in the same line, to allow the output of the first command to be used as the input to the second.
```bash
# List only the first 20 lines of only samples SRR1039508 (col 2) and SRR1039523 (col 17)
cut -f 1,2,17 all_counts.txt | head -n 20
```

Similarly to how we used the pipe operator (|) above, we could use the redirect operator(>) to send the output of the cut command to create a new counts file, that only contains the columns 1 (gene IDs), and samples in columns 2 and 17.
```bash
# Print the counts from SRR1039508 and SRR1039523 to a new file
cut -f 1,2,17 all_counts.txt > all_counts_sub.txt

# look at head of this new file
head all_counts_sub.txt
```

### Pattern matching with *grep*

Often we will want to pull a specific piece of information from a large file, let's say that we were interested in the read counts for a specific gene, ALDH3B1 (Ensembl ID: ENSG00000006534). We can use the `grep` command to search for this ID, or any other character string we are interested in, in our counts matrix.
```bash
# Get the count data for ENSG00000006534 (ALDH3B1) from all_counts.txt
grep "ENSG00000006534" all_counts.txt
```

`grep` is a pattern recognition tool that searches in files for a character string we can define. We can define the entire character string, as we did above, or combine regular characters with special characters (or 'wildcards') to search for specific types of matches. Some commonly used special characters are included in the table below.

Operator | Effect
---|---
\* | wildcard stands for any number of anything
^ | start of the line
$ | end of the line
[0-9] or \d| any number (0123456789)
[a-z]| any lowercase letter
[A-Z]| any uppercase letter
\t | a tab


These regular expressions can be used with any of the tools that you have learned thus far, so if we wanted to list all of the files in our directory that end in .txt we could use the following command.

```bash
# List all files that end in .txt
ls *.txt
```

We can even enhance the power of these regular expressions by specifying how many times we expect to see the regular expression with quantifiers.

Quantifier| Operation
---|---
X* | 0 or more repetitions of X
X+ | 1 or more repetitions of X
X? | 0 or 1 instances of X

Now let's use some of these regular expressions in a `grep` command  to see their utility. Let's use regular expressions to see how many genes have zero reads counted for the first four samples. The flag `-P` indicates that we will be using perl-style regular expressions in the pattern we are searching for, you can use `grep --h` to learn more about available flags for the `grep` command.

```bash
# Count the number of genes with no reads in the first four samples
grep -P "^ENSG[0-9]*\t0\t0\t0\t0\t" all_counts.txt| wc -l

# Count the number of genes with no reads expressed in any of the samples
grep -P "^ENSG[0-9]*\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0$" all_counts.txt| wc -l
```
> As we move through the subsequent lessons, we will introduce more complex bash commands in order to manipulate common bioinformatics file types. If you are ever confused about what a command does, remember you can always use `man` to check out the manual page (or Google it). It you are confused about how commands are used in conjunction with each other, it can also be helpful to break them down and run parts individually, in order to understand what the constituent parts do.

### Further learning

- Check out the Bash/Unix cheat sheet links in the GitHub directory, and try out a few other commands on the `all_counts.txt` file. Use the `man` command to learn more about these commands and how they should be used.
