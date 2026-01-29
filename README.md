This repository contains additional material to help with the usage of
[notame](https://bioconductor.org/packages//release/bioc/html/notame.html),
[notameStats](https://bioconductor.org/packages//release/bioc/html/notameStats.html),
and
[notameViz](https://bioconductor.org/packages//release/bioc/html/notameViz.html).
Currently, it contains the following files:

- `notame_workflow.R` - A workflow script for the usage of notame.
  (preprocessing), notameStats (statistics), and notameViz (visualisations) with
  example data.

## How to Run the Script?

1. Click on the green button "Code" in the repository and download the zip
   archive.
2. Unzip the archive.
3. [Install](https://cloud.r-project.org/) R programming language on your
   system.
4. [Install](https://posit.co/downloads/) RStudio.
5. Open `notame_workflow.R` in RStudio. Ensure that you open the file in the
   folder where you unzipped the archive (the folder should contain the
   `renv.lock` file).
6. Run the lines up until and optionally including `renv::restore()`. This step
   should install all the required packages to be able to run the script.
