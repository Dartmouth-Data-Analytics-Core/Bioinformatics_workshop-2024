# Welcome to the DAC Fundamentals of Bioinformatics workshop 

Before you attend the workshop there are a couple of things we would like you to do to get setup for using the tools that will be required during the course of the workshop. Please read through each of the sections below to ensure you are prepared to attend the workshop. We strongly recommend that you run through these steps several days in advance of the workshop, in case any troubleshooting is required.


## The terminal emulator
---

A terminal emulator is a more streamlined way of navigating a computational environment (once you get the hang of it). We will cover some basic commands to help orient you with using the terminal to interact with your machine on Day 1 of the workshop.

If you are using a Mac there is a terminal emulator program pre-installed on your machine, if you are using another OS we have included some links to popular terminal emulator programs below. Please select one of them download the program and open it up for the next step.

Operating system| Terminal emulators
---|---
Mac| Terminal (comes pre-installed)
Windows| [MobaXterm](https://mobaxterm.mobatek.net/download.html) <br> [PuTTY](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html)
Linux| Konsole, Terminal, etc. (should be pre-installed but depends on the desktop environment you are running)

## VPN client 
---

For those that are attending the workshop off campus you will need to download the VPN client and log in to access the discovery cluster, which we will use for the lessons on Days 1 and 2. You can download the VPN client [here](https://services.dartmouth.edu/TDClient/1806/Portal/KB/ArticleDet?ID=66806) and you will be able to log onto the VP network with your Dartmouth credentials. If you've registered for the workshop but do not have Dartmouth credentials please reach out to Shannon.Soucy@dartmouth.edu to get a sponsored account. 


## The discovery HPC system ##
---

For those of you that indicated that you did not have an account on *discovery* you will need to request one [here](https://rc.dartmouth.edu) **AT LEAST** two (business) days before the workshop begins. There is a green button at the top of the page that says **Request an Account**, once you click the button you will be prompted to log in with your netID and password and then you can fill in the form to request your account.

Once you have a discovery account you can follow along with this video [here](https://youtu.be/VoHBlblsQfg) to log onto discovery using the command line.

To log onto discovery we will use the secure shell command `ssh`. 

```bash
ssh netID@discovery.dartmouth.edu
```

You will be prompted for a password and when you start typing nothing will show up in the prompt, I assure you though your keystrokes are recorded and you will be logged onto the discovery HPC environment if your password is correct. If your password is incorrect you will receive the warning "Permission denied, please try again." and will be prompted to enter your password again.



## Installing an SFTP client ##
---

For those of you that are new to the command line this might be an easier way to move files between the HPC environment and your local machine. An SFTP client stands for secure file transfer protocol and will enable you to drag and drop files as you might in a finder window between your local machine and a remote location. 

We recommend FileZilla, which works on Mac, Windows, and linux operating systems. You can download [FileZilla](https://filezilla-project.org/download.php?show_all=1) by following the link and selecting the version that is correct for your OS, then open the program to ensure that you have downloaded it successfully. Once you have Filezilla installed you can use this [video](https://youtu.be/A8w8Uw1OILA) to guide you in linking the SFTP client to your account on discovery.

<img src="figures/filezilla.png" height="80%" width="80%" />


## Install the Integrative Genomics Viewer (IGV)
---

We will be using the [Integrative Genomics Viewer (IGV)](http://software.broadinstitute.org/software/igv/), a genome browser produced by researchers at the Broad Institute, to explore and visualize genomics data on Day 2. 

<img src="figures/igv.png" height="100" width="100"/>

You will need to download and install the IGV desktop application for your operating system before the workshop begins. The latest versions of IGV can be found at their [downloads page](http://software.broadinstitute.org/software/igv/download). After installing IGV, try opening the application on your computer to confirm the installation was successful. Please ensure that the **Human (hg38)** genome is loaded in the IGV version you downloaded.


## Setting up an R project ##
---

We will be using R-Studio to explore and analyze genomics data on day 3, therefore we ask that you have R and R-Studio installed prior to attending the workshop. You will need to be running at least version 3.6 to ensure all of the packages needed will run smoothly. The latest versions for R and R-Studio can be found [here](https://cran.r-project.org) and [here](https://rstudio.com/products/rstudio/download/).

Next you will need to set up a new project in R-Studio for this workshop. Projects in R are like containers for various jobs that you will perform, the history of the project will be loaded when you open a new project. By using a project you can install all of the tools ahead of time and they will be there for you when you need them during the workshop. In R-Studio under File select New directory and then select New Project and name the project something you will remember (bioinfo_workshop).

Now that you have your project loaded, run the following code to install all of the packages we will be using during the workshop. For those of you that haven't used RStudio before we have made a [video](https://youtu.be/UtZHS-q7buI) showing the successful installation of the R packages you will need using the commands below. 

I bumbled the description of the code chunks with the nested loops in the video, so here is a better description for those that are interested: 

There are two loops and each loop starts with an if statement. The first loop states "if `biomaRt` is not installed enter this loop" and the second one "if `BioCManager` is not installed enter this loop", when the condition is fulfilled (the package is *not* installed) the loop is entered and the function `install.packages` is used to install the package. Each loop is exited once the packages are installed and the package is loaded with the 'library' function to make the functions contained in the package available during the current session.


```r
if (!any(rownames(installed.packages()) == "biomaRt")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("biomaRt")
}
library(biomaRt)

if (!any(rownames(installed.packages()) == "IRanges")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("IRanges")
}
library(IRanges)

if (!any(rownames(installed.packages()) == "GenomicRanges")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("GenomicRanges")
}
library(GenomicRanges)

if (!any(rownames(installed.packages()) == "Gviz")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("Gviz")
}
library(Gviz)

if (!any(rownames(installed.packages()) == "org.Hs.eg.db")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)

if (!any(rownames(installed.packages()) == "EnsDb.Hsapiens.v86")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("EnsDb.Hsapiens.v86")
}
library(EnsDb.Hsapiens.v86)

if (!any(rownames(installed.packages()) == "GenomicFeatures")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("GenomicFeatures")
}
library(GenomicFeatures)

if (!any(rownames(installed.packages()) == "VariantAnnotation")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("VariantAnnotation")
}
library(VariantAnnotation)

if (!any(rownames(installed.packages()) == "TxDb.Hsapiens.UCSC.hg38.knownGene")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
}
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

if (!any(rownames(installed.packages()) == "TxDb.Mmusculus.UCSC.mm10.knownGene")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
}
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

if (!any(rownames(installed.packages()) == "BSgenome")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("BSgenome")
}
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

if (!any(rownames(installed.packages()) == "ChIPseeker")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("ChIPseeker")
}
library(ChIPseeker)

BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")

BiocManager::install("BSgenome.Mmusculus.UCSC.mm10.masked")

sessionInfo()
```


When installing R packages using the code chunk above a lot of red text prints to the screen, most of the text are messages about the packages and dependencies being installed. However some of these messages could be errors indicating a package was not installed. It is good practice to read through these messages however a faster way to identify if all the software you need is installed is to run the following code. This code first creates a list of the packages you need called `pkg.list`, next we loop through the list (more on this on Day 2) to check if each package was successfully installed. If the package is NOT installed the name of the package will print to your screen and you can find the install command in the code chunk above and re-run that code chunk. 

```
pkg.list<-c("biomaRt","IRanges","GenomicRanges","Gviz","org.Hs.eg.db","EnsDb.Hsapiens.v86","GenomicFeatures","VariantAnnotation","TxDb.Hsapiens.UCSC.hg38.knownGene","TxDb.Mmusculus.UCSC.mm10.knownGene","BSgenome","ChIPseeker","BSgenome.Mmusculus.UCSC.mm10","BSgenome.Mmusculus.UCSC.mm10.masked")

# loop to check which packages are NOT installed 
for(i in 1:length(pkg.list)) {
  if (!any(rownames(installed.packages()) == pkg.list[i])){
    print (pkg.list[i])
  }
}
# If a package name is printed to the screen you must install that package before continuing with the workshop. 
```


## Downloading the data 
---

The commands that you will be following can be found in markdown `(.md)` files where there is a brief description of each command and how it is applied to the data and what it does followed by an example command that you can copy and paste into the terminal window. Day 1 and 2 will be using the terminal window on your local machine, with an open `ssh` connection to discovery, as we will be running `bash` code. For some of day 2 and most of day 3 you will be using RStudio on your local machine to run the commands in the markdown files (`.md`) located in this GitHub repo. 


In your terminal window **on your local machine** navigate to where you want to download the files needed for this workshop. 

**On Monday** before you log onto the first zoom session we will make the workshop materials public and you should download those to your local machine (preferably in the same location as you downloaded the setup materials) with the following command: 

```bash
git clone https://github.com/Dartmouth-Data-Analytics-Core/Bioinformatics_workshop-Dec-2022/
```
If you run this command before Monday morning you may not download the most up to date lessons and your materials may not be identical to those that are presented in the workshop.

---

If you have issues with any part of the installation and setup, please reach out to us directly DataAnalyticsCore@groups.dartmouth.edu, or attend optional bioinformatics help hours **Friday, December 2, 2021 12:30-1:30PM** the zoom link can be found in your welcome email. 

