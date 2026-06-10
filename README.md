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

### Docker Container (For Advanced Users)

This repository also includes a Docker container with all the packages
pre-installed. Learn how to install it on
[Linux](https://docs.docker.com/engine/install/).

Install [rootless Docker](https://docs.docker.com/engine/security/rootless/).

After that, in the terminal run

```shell
docker run --rm -it -e PASSWORD=1234 -p 8787:8787 \
  ghcr.io/hanhineva-lab/notametutorials:latest
```

In your browser, go to `http://localhost:8787` and enter RStudio with the
username `root` and the password `1234`. In the R console, run
`setwd("/project")` to enter the working directory of the script.

> [!TIP]
> If you modify anything inside the Docker container and then stop it, the data
> will be lost. To persist it, add the `--volume ./:/project` parameter to map
> the current directory on your computer to the `./project` directory inside
> the Docker container. Read more about the volumes in the [official
> documentation](https://docs.docker.com/engine/storage/volumes/).
