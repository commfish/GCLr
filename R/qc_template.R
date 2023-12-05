#' Quality Control Analysis
#' 
#' This function opens the `qc_template.rmd` file. The details section of this document walks the user through the process of running the lab quality control analysis script.
#' 
#' @details If you are just running the project data analysis, walk through steps 1-14 in this guide.\cr If you are also running the QC data analysis, walk though steps 1-18 in the guide.
#' 
#' 1) Open R Studio (If you are reading this, you probably RStudio open already `r emo::ji("silly")`) 
#' 2) - **Option 1:** From the menu bar, select: `File -> New File -> R Markdown` and then select “From Template” on the left menu; “qc_template” should be highlighted. Hit “OK” \cr 
#'           \figure{NewRMDWindow_small.PNG}
#'    - **Option 2:** Simply run `qc_template()` in your console (see tab marked "Console" in RStudio) and the file will open (i.e., paste `qc_template()` in your console and hit return)
#' 3) At the top of the script, fill in the following: 
#'    - `title:` “Project number”
#'    - `author:` “Your name”
#'    - `date:` The date auto fills
#' 4) “Save as” the file, and it will pop up in your QC folder (if it doesn’t, check your working directory). Close RStudio. Double-click on the R file that just showed up in the QC folder and work in there. (This is an important step because it sets your working directory to the folder where you saved the script.)
#' 5) Click the green arrow \figure{RunChunkGreenArrow.png}(aka "Run Current Chunk" button) to run the first code chunk before the Intro. Code chunks are the gray boxes containing code. 
#' 6) Scroll down to the `input` code chunk and fill in the following:
#'    - `dirqc <- “paste the project directory for the QC folder”` the directory path must have 2 backslashes (`\\`) or 1 forward slash (`/`) between folders or it won't work. A simple way to get this is to run `getwd()` in the console, then copy and paste the directory path from the console into the script.
#'    - `species <- “species”` it must be lowercase.
#'    - `project <- “Project number”`
#'    - `project_type <- “Just copy paste one of the platforms”`
#' 7) Click \figure{RunChunkGreenArrow.png} in the `inputs` code chunk.
#' 8) Click \figure{RunChunkGreenArrow.png} in `optional_inputs` code chunk.
#' 9) In the Notebook prep section, click \figure{RunChunkGreenArrow.png} in the `knitprep` code chunk, then
#' 10) Click \figure{RunChunkGreenArrow.png} in `formatfun` code chunk.
#' 11) Scroll down to the `Read Project Genotypes` section and click \figure{RunChunkGreenArrow.png} in the `project_genos` code chunk. This could take a few minutes if it is a large project; all genotypes are being pulled.
#' 12) Keep moving down the script and click on \figure{RunChunkGreenArrow.png} in each of the remaining code chunks. You must click on all the green arrows in order as they build upon each other.
#'      Here are the remaining code chunks and a brief description of what they do:
#'      - `fail_rate` - calculate the project failure rate
#'      - `fail_rate_silly_table` - calculate failure rate by silly (aka collection)
#'      - `fail_rate_silly_plot` - plot the failure rate by silly
#'      - `fail_rate_locus` - plot the failure rate by locus
#'      - `fail_rate_plate_table` - create a table showing the failure rate by plate for the project
#'      - `fail_rate_plate_plot` - plot the failure rate by plate
#'      - `fail_rate_overall` - calculate the overall failure rate for the project
#'      - `samp_size_locus` - get the number of samples by locus for each silly
#'      - `samp_size` - get the number of fish genotyped for the project
#'      - `alt_spp` - check for alternate species (i.e., if the wrong species was extracted). Note that this is only done for chum and sockeye and if there are certain markers in the markersuite. Otherwise, NA is added to the sample size table under `wrong_spp`.
#'      - `samp_size_alts` - take a sneak peek at the sample sizes table after the alternate species samples size is added 
#'      - `alt_spp_table` - view the table containing failed and alternate markers for each fish.
#'      - `miss_loci` - check for project individuals missing \>=20% of genotypes (i.e. the 80% rule).
#'      - `samp_size_miss` - check out the sample size table, after removing the fish with poor quality data
#'      - `dupcheck` - Remove duplicate individuals within the same collection. Typically we specify duplicates as a pair of individuals that share >=95% of genotypes. 
#'      - `samp_size_dups` - check out the sample size table, after checking for duplicates
#'      - `samp_size_final` - this table shows the final sample size for each silly 
#'      - `select_qc_fish` - select high-quality fish for the QC analysis. A popup window will ask you if you want to select QC fish yes/no **Note: there will be a function to do this in the future, but it's not ready for prime time, so disreguard this code chunk.** 
#'        - *This ends the project data analysis portion of the script.*
#'  13) Manually enter a project data analysis summary:
#'             - `Name:` Your name
#'             - `Date:` Enter the date
#'             - `Comments:` Anything you noted from project review.
#'             
#'        **At this point, you can simply save the file and exit or continue onto the QC data analysis section.**
#' 
#'  14) If you aren't continuing with the QC analysis, click save then close. When closing, a popup window will ask you if you want to save, click "Don't Save".
#'        **If you are moving on to the QC data analysis follow steps below.**
#'        
#'  15) Run each code chunk (click \figure{RunChunkGreenArrow.png}) in the *QC Data Analysis section:
#'      - `qc_genos` - read in the QC genotypes
#'      - `concordance_files` - read in the concordance files
#'      - `conflict_plate` - take a look at the number of conflicts by plate
#'      - `conflict_silly` - take a look at the number of conflicts by silly
#'      - `conflict_locus` - take a look at the number of conflicts by locus
#'      - `qc_reruns` - calculate sample size by locus for QC genotypes and look a list of failed sillys (i.e., less than 80% successful genotypes)
#'      - `plot_successrate` - create an interactive plot of success rates
#'      - `qc_rm_miss` - check for QC individuals missing \>=20% of genotypes (i.e. the 80% rule)
#'      - `qc_conflict_rate` - get QC individual conflict rate and plot. If there are no conflicts, no plot will be produced.
#'      - `qc_dupcheck` - identify anything to investigate by running a duplicate check on anything with high degree of conflicts
#'      - `dup_check_view` - take a look at the duplicate check results
#'      - `summ_by_silly` - table overall summary by silly for the project (i.e., failure rates, sample sizes)
#'      - `conf_by_silly` - table conflicts by silly for the project
#'      - `conf_by_locus` - table conflicts by locus for the project
#'      
#'  16) Manually enter a QC data analysis summary:
#'             - `Name:` Your name
#'             - `Date:` Enter the date
#'             - `Comments:` Anything you noted from project review.
#'      
#'      
#'  17) Run `savedata` code chunk to save workspace image
#'  18) Save file: `File -> Save`. Save the file in the QC directory and rename the notebook html file to something informative (e.g.,  *project#*_QC_final.nb.html)
#'      
#' @note Keep track of what tables and figures you find informative and note the ones that do not bring you joy so we can remove/hide them in the next package version.
#' 
#' @export
qc_template <- function(){
  
  template <- system.file("rmarkdown/templates/qc_template/skeleton", "skeleton.rmd", package = "GCLr") 
  
  file.copy(from = template, to = path.expand("~/qc_template.rmd"))
  
  shell(path.expand("~/qc_template.rmd"))
 
}