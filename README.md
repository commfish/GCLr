
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GCLr <img src="man/figures/logo.png" align="right" height="100" /></a>

This is a temporary README for the GCLr package. This document is
intended as a guide help collaborators add package documentation to
their assigned functions. The final README will be added after all
functions have been documented.

## Introduction

The purpose of this guide is to provide a brief tutorial on how to
modify existing functions for inclusion in an R package. Below we
highlight some of the commonly used syntax used for function
documentation and provide steps to assist with documentation. While
functions can be modified manually, we’re recommending ChatGPT. ChatGPT
seems to work well and will help with standardizing the language, uses,
etc. Further, it knows about additional formatting options that make the
documentation look nice. Also, it will save you time! In short, you’ll
prompt ChatGPT to re-write the header of the function and it will do so.
Notice that I’m only asking for the header/documentation. It will
re-write entire functions, but 1) we don’t want this without a thorough
review of the new function, 2) there is a character limit, so you may
not get the whole thing, and 3) it may make things up. Note that ChatGPT
answers require manual intervention, but this is what we want. We’re
using it as a tool to get the barebones documentation, before pushing
the files back to the repo.

### Repository

The new repo is here: <https://github.com/commfish/GCLr> 

And functions are located here: <https://github.com/commfish/GCLr/tree/main/R>

These are the latest and greatest files and what you should be working
on (i.e., not Kyle’s GH). They have been renamed and a few lines were
removed, where necessary. They were uploaded to serve as a starting
point, so we have records of the changes.

To edit the package functions, you'll have to clone the repository to your local computer. 
Assuming you have rstudio and git, this can be done in RStudio via:

    File > New Project > Version Control > Git
        repository URL: https://github.com/commfish/GCLr/
        Project directory name: (this should autofill but you can call it whatever you want)
        Create project as subdirectory of: (where you want the files to go)_

You should now see the project and can edit the files.

### Installation

You can install the package using devtools. After editing files, and
pushing back to GitHub, it’s recommended to install again so you can see
if your function broke things. Before you reinstall, it’s best to shut
down R, delete the package folder in your library folder, and finally
reinstall it all.

``` r

install.packages("devtools")

library(devtools)

devtools::install_github("commfish/GCLr")
```


## Roxygen syntax

In R package documentation, there are various tags that can be used to
provide additional information about the function or package. They are
identified via tags, and you can also use standard markdown formatting
to make the documentation look cleaner, easier to read, etc.

### Tags

This section describes some of the commonly used tags, a few of which
we’ve marked as mandatory. Feel free to include some or all of them
where applicable/known but we’re hoping to include all the “mandatory”
tags for each function. Note that there are more tags, so feel free to
explore and let everyone know if you find something useful! Here is a
link to the some of the tags:

1.  @title: Mandatory; Specifies the title or name of the function. This
    is the first line of the file. Note use of “@title” is optional

2.  @description: Mandatory; Provides a brief description of the
    function or package. This the second line(s) of the file. Note use
    of “@description” is optional

3.  @param or @arg: Mandatory; Describes function parameters or
    arguments, including their names, types, and descriptions.

    1.  ChatGPT does decent here but check to make sure they are right,
        and all accounted for
    2.  Each input goes on a separate @param line

4.  @return: Mandatory; Describes the return value or values of the
    function.

5.  @export: Mandatory; Indicates that the function should be exported
    and made available to users of the package.

    1.  Make sure this is here or package build will fail. Sometimes
        ChatGPT misses this.
    2.  Add it to the end of the header.

6.  @details: Mandatory if old function has warnings, notes, etc.;
    Provides detailed information about the function or package.

7.  @family: Assigns the function to a specific family or group. This
    will be added later.

    1.  You can e decided to assign family to applicable functions so
        that they are easier to find, search for, etc. This is kind of
        like our original rationale for multiple packages.
    2.  You can have multiple families for each function.

8.  @examples: Provides examples of how to use the function.

    1.  Make sure this looks okay (e.g., spacing, language, etc.), but
        no need to spend a ton of time here since we will be adding set
        examples throughout.
    2.  Only toss in old example, do not modify further. This will be
        done later.

9.  @seealso: Specifies related functions or topics that users may want
    to explore. Note that this is generally for internal links (i.e.,
    other functions in the package)

10. @references: Includes references to articles, books, or other
    resources related to the function or package. Note that this is
    generally external links (i.e., not other functions in the package)

11. @keywords: Specifies keywords or tags associated with the function.

12. @import: Remove if added by ChatGPT; Specifies packages that need to
    be imported for the function to work.

13. @importFrom: Remove if added by ChatGPT; Specifies specific
    functions from other packages that need to be imported.

14. @Author: Specifies who wrote the function.

### Tips

The purpose of this section is to outline a few additional tricks that
we’ve found helpful in documenting the functions. Note that some of
these are stylistic (i.e., optional) and others should be done, or the
package may not build properly.

1.  Mandatory; Be careful with the two lines at top. ChatGPT sometimes
    doesn’t quite get these. It may depend on existing spacing or line
    breaks in the header.
    1.  The first line is a very brief (i.e., handful of words)
        description of the function (i.e., @title, in title case)
    2.  The second is a more in-depth description, where a sentence or
        two is okay (i.e., @description).
2.  Mandatory; Remove any instances of loading packages or requiring
    packages in the beginning.
    1.  These should not be included since the package takes care of
        this.
    2.  I already removed the “pacman” lines from all functions but
        there may be remnants or other cases (e.g.,
        while(!require(tidyverse)){ install.packages(“tidyverse”)})
3.  Mandatory; Remove any @import or @importFrom tags – these are useful
    but we’re going a different route with the package and keeping them
    tosses a bunch of warnings.
4.  Optional; Use square brackets when referring to functions in the
    documentation: \[GCL::test_LD()\]
    1.  This will link relevant functions and provide more information
        when referencing in the documentation.
    2.  An example can be seen in – \[GCL::summarize_LD()\]
5.  Optional; Use "_\itemize{_" to generate bulleted lists
    1.  These make the documentation cleaner and easier to follow
    2.  An example can be seen in – \[GCL::confusion_matrix()\] - with a truncated version here:

``` r

#' @return Named list of tibble(s) in long format:
#'   \itemize{
#'     \item \code{group_group}: A tibble containing the following variables:
#'       \itemize{
#'         \item \code{repunit}: The known repunit (aka reporting group).
#'         \item \code{inferred_repunit}: The inferred repunit.
#'         \item \code{mean_group_group_scaled_like}: The average probability of each individual from a group originating from a pop in that group.
#'       }
#'     \item \code{pop_group}: A tibble containing the following variables:
#'       \itemize{
#'         \item \code{collection}: The known population.
#'         \item \code{inferred_repunit}: The inferred repunit.
#'         \item \code{mean_pop_group_scaled_like}: The average probability of each individual from a pop originating from a pop in that group.
#'       }
#'  } 
```

## Example Using ChatGPT

Now, we’ll demonstrate an example workflow using ChatGPT. If you don’t
have one, you can sign up for a free account here
<https://chat.openai.com/>. Then you can start a new chat and ask to
modify the function. The best way to do this is to copy the entire
function and paste it into the chat.

Here is a suggested prompt – 
***Convert the function documentation so it can be included in an R package. “paste the function here”***

Here are the general steps used to convert add_tree_color.r

1)  Open the function locally or on the package repo

2)  Select all via Control + A

3)  Copy via Control + C

4)  Switch to ChatGPT

5)  Type in the prompt above and paste the function into the chat window
    via Control + V

6)  Hit enter.

7)  Watch it go! a. You can click “stop generating”, if it continues
    past the header, since we don’t care about the function.

8)  Add the new header into the file and manually edit it, fixing any
    issues.

9)  Once finished, save.

Here is the code from the first attempt – this looks reasonable but
there are still things to fix. For example, @example is not quite right,
some @params have funny spaces, the description is spaced funny, etc.

``` r

#' Add Group Colors to Phylogenetic Tree and Modify Tip Labels
#'
#' This function takes a tree object produced by `ape::nj()` and adds group colors
#' to the branch lengths for each reporting group. Additionally, it can be used to
#' modify the tip labels with new names and optional symbols.
#'
#' @param tree A phylogenetic tree object (class = "phylo") produced by the ape package.
#' @param currentnames A character vector of the tip labels (pop names) in the tree object.
#' @param treenames A character vector of new tip labels (pop names) to replace current names.
#' @param groupvec A numeric vector indicating the group affiliation for each population.
#'   Note that `currentnames`, `treenames`, and `groupvec` should be of the same length and order.
#' @param regioncol A vector of colors with the same length as `max(groupvec)`,
#'   representing the number of groups. Colors can be specified as hexadecimal,
#'   R color names, or the number corresponding to an R color name in `colors()`.
#' @param regionpch A vector with the same structure as `regioncol`, where each element
#'   is an R pch (plot character) number.
#' @param write_nexus Logical indicating whether to write a Nexus tree file that can be opened in FigTree.
#'   Default is `FALSE`.
#' @param file If `write_nexus = TRUE`, the full file path for writing out a Nexus tree file,
#'   including the .nex extension.
#'
#' @return A new tree object with colors and modified tip labels.
#'
#' @examples
#' # Example usage
#' tree <- ape::nj(dist.matrix)
#' currentnames <- tree$tip.label
#' treenames <- c("newname1", "newname2", "newname3") # Replace with actual new names
#' groupvec <- c(1, 1, 2) # Replace with actual group affiliations
#' regioncol <- c("red", "blue") # Replace with actual colors
#' colortree <- add_tree_color(tree, currentnames, treenames, groupvec, regioncol)
#'
#' @references
#' - FigTree: Tree Visualization and Annotation: https://github.com/rambaut/figtree/releases
#' - R color names: See `colors()` function in R.
#'
#' @export
add_tree_color <- function(tree, currentnames, treenames, groupvec, regioncol, regionpch = NULL, write_nexus = FALSE, file = NULL) {
  # Function code goes here
}
```

Once you’ve completed the function, push it back to the repository and
check it off your list on the repo issue. Please add comments or
anything else you think would be useful to your tasks.

# Function Check
To view the documentation of the newly generated functions, and ensure it's working, you can use:

``` r
devtools::document()
```

## Big picture steps

1)  Copy function

2)  Paste into ChatGPT

3)  Copy documentation from ChatGPT output

4)  Paste it above the function in R

5)  Edit to fit the guidelines described here (removing old
    documentation)

6)  Save file

7)  Do a fresh pull from GitHub

8)  Commit with meaningful message (or just repeat 1 – 7 then carry on
    with 8 – 10 for all of them)

9)  Mark task completed (check the box) in your issue

10) Move to the next function
