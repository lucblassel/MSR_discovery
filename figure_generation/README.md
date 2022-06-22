Mapping-friendly sequence reductions: going beyond homopolymer
compression
================
Luc Blassel, Paul Medvedev & Rayan Chikhi
Jun 22 2022

-   [Setup](#setup)
    -   [Needed functions](#needed-functions)
    -   [Data loading and
        pre-processing](#data-loading-and-pre-processing)
-   [Main text](#main-text)
    -   [Figure 1](#figure-1)
    -   [Figure 2](#figure-2)
    -   [Figure 3](#figure-3)
    -   [Figure 4](#figure-4)
    -   [Figure 5](#figure-5)
    -   [Table 1](#table-1)
-   [Supplementary Material](#supplementary-material)
    -   [Additional functions for
        sup-mat](#additional-functions-for-sup-mat)
    -   [Table B.1](#table-b1)
    -   [Figure C.1](#figure-c1)
    -   [Figures D.1 to D.5](#figures-d1-to-d5)
    -   [Figure E.1](#figure-e1)

This notebook was used to generate the table and figures used in the MSR
paper *([DOI](placeholder))*.

# Setup

Necessary libraries and general settings

``` r
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(kableExtra)
library(pafr)
library(rjson)
library(ggridges)

# Setting general theme options 
ggplot2::theme_set(ggplot2::theme_classic())
ggplot2::theme_update(
  strip.text=element_text(size=12, hjust=0, vjust=1),
  axis.title=element_text(size=12, hjust=0.5, vjust=0.5),
  axis.text=element_text(size=11, hjust=0.5, vjust=0.5),
)

# Mapq thresholds to show as points
point_thresholds <- seq(0, 60, 10)

# Undesirable MSRs 
to_reject <- c('R(3,False,60,1)',
 'R(3,False,5,2)',
 'R(3,False,4,2)',
 'R(3,False,48,2)',
 'R(3,False,38,0)',
 'R(4,True,27,6)',
 'R(4,True,24,10)')

colour_levels <- c("HPC", "raw", "MSR", "MSR E", "MSR F", "MSR P")
colors <- c("#F8766D", "#00BA38", "#619CFF26", "#004c6d", "#0098bb", "#51e3ff")
to_show <- c("MSR_E", "MSR_F", "MSR_P")
```

## Needed functions

### Custom Facet-wrapping

``` r
##########################
# SCALE OVERRIDING UTILS #
##########################

# from: https://dewey.dunnington.ca/post/2018/modifying-facet-scales-in-ggplot2/

scale_override <- function(which, scale) {
  if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
    stop("which must be an integer of length 1")
  }
  
  if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
    stop("scale must be an x or y position scale")
  }
  
  structure(list(which = which, scale = scale), class = "scale_override")
}

CustomFacetWrap <- ggproto(
  "CustomFacetWrap", FacetWrap,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    # make the initial x, y scales list
    scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)
    
    if(is.null(params$scale_overrides)) return(scales)
    
    max_scale_x <- length(scales$x)
    max_scale_y <- length(scales$y)
    
    # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
    for(scale_override in params$scale_overrides) {
      which <- scale_override$which
      scale <- scale_override$scale
      
      if("x" %in% scale$aesthetics) {
        if(!is.null(scales$x)) {
          if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
          scales$x[[which]] <- scale$clone()
        }
      } else if("y" %in% scale$aesthetics) {
        if(!is.null(scales$y)) {
          if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
          scales$y[[which]] <- scale$clone()
        }
      } else {
        stop("Invalid scale")
      }
    }
    
    # return scales
    scales
  }
)

facet_wrap_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  # sanitize scale overrides
  if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) || 
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }
  
  facet_super$params$scale_overrides <- scale_overrides
  
  ggproto(NULL, CustomFacetWrap,
    shrink = facet_super$shrink,
    params = facet_super$params
  )
}
```

### Table formatting functions

``` r
fmt_perc <- function(num) {
  if (is.na(num)) {
    return ("NAN")
  }
  if (num >= 0) {
    return(
      paste0("+", format(num, digits=1, width=3, justify="right"), "%")
    )
  }
  return (
    paste0(format(num, digits=1, width=3, justify="right"), "%")
  )
}

fmt_vec <- function(vec) {
  return(
    lapply(vec, fmt_perc)
  )  
}
```

## Data loading and pre-processing

``` r
evals <- read.csv("./data/MSR-evals.csv") %>%
    filter(!(renamed %in% to_reject)) %>%
    mutate(col=str_replace(MSR, "_", " ")) %>%
    mutate(col=factor(col, levels=colour_levels))
```

# Main text

## Figure 1

|     ![Figure 1](figures/figure1.svg)      |
|:-----------------------------------------:|
| Figure 1 was made by hand using inkscape. |

## Figure 2

|     ![Figure 2](figures/figure2.svg)      |
|:-----------------------------------------:|
| Figure 2 was made by hand using inkscape. |

## Figure 3

|                           ![Figure 3](figures/figure3.svg)                            |
|:-------------------------------------------------------------------------------------:|
| Figure 3 is a modified version of panel A of figure 5. It was modified with inkscape. |

## Figure 4

|     ![Figure 4](figures/figure4.svg)      |
|:-----------------------------------------:|
| Figure 4 was made by hand using inkscape. |

## Figure 5

<img src="figures_files/figure-gfm/plot_custom_reordered, , -1.png" style="display: block; margin: auto;" />

## Table 1

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="2">
</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="4">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Whole HG

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="4">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Whole HG

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="4">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Droso

</div>

</th>
</tr>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="2">
</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="4">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

minimap2

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="4">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

winnowmap

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="4">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

minimap2

</div>

</th>
</tr>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="2">
</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Frac

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Err

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Frac

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Err

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Frac

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

Err

</div>

</th>
</tr>
<tr>
<th style="text-align:left;">
MSR
</th>
<th style="text-align:right;">
t
</th>
<th style="text-align:left;">
fWM
</th>
<th style="text-align:left;">
PfWM
</th>
<th style="text-align:left;">
eWM
</th>
<th style="text-align:left;">
PeWM
</th>
<th style="text-align:left;">
fWW
</th>
<th style="text-align:left;">
PfWW
</th>
<th style="text-align:left;">
eWW
</th>
<th style="text-align:left;">
PeWW
</th>
<th style="text-align:left;">
fDM
</th>
<th style="text-align:left;">
PfDM
</th>
<th style="text-align:left;">
eDM
</th>
<th style="text-align:left;">
PeDM
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
HPC
</td>
<td style="text-align:right;">
60
</td>
<td style="text-align:left;">
0.935
</td>
<td style="text-align:left;">

-   0%
    </td>
    <td style="text-align:left;">
    1.85e-03
    </td>
    <td style="text-align:left;">

    -   0%
        </td>
        <td style="text-align:left;">
        0.894
        </td>
        <td style="text-align:left;">

        -   0%
            </td>
            <td style="text-align:left;">
            1.43e-03
            </td>
            <td style="text-align:left;">

            -   0%
                </td>
                <td style="text-align:left;">
                0.957
                </td>
                <td style="text-align:left;">

                -   0%
                    </td>
                    <td style="text-align:left;">
                    2.27e-03
                    </td>
                    <td style="text-align:left;">

                    -   0%
                        </td>
                        </tr>
                        <tr>
                        <td style="text-align:left;">
                        raw
                        </td>
                        <td style="text-align:right;">
                        60
                        </td>
                        <td style="text-align:left;">
                        0.921
                        </td>
                        <td style="text-align:left;">
                        -1%
                        </td>
                        <td style="text-align:left;">
                        1.86e-03
                        </td>
                        <td style="text-align:left;">
                        +0.2%
                        </td>
                        <td style="text-align:left;">
                        0.932
                        </td>
                        <td style="text-align:left;">

                        -   4%
                            </td>
                            <td style="text-align:left;">
                            1.75e-03
                            </td>
                            <td style="text-align:left;">

                            -   23%
                                </td>
                                <td style="text-align:left;">
                                0.958
                                </td>
                                <td style="text-align:left;">
                                +0.04%
                                </td>
                                <td style="text-align:left;">
                                2.27e-03
                                </td>
                                <td style="text-align:left;">
                                -0.04%
                                </td>
                                </tr>
                                <tr>
                                <td style="text-align:left;">
                                MSR_F
                                </td>
                                <td style="text-align:right;">
                                50
                                </td>
                                <td style="text-align:left;">
                                0.938
                                </td>
                                <td style="text-align:left;">
                                +0.4%
                                </td>
                                <td style="text-align:left;">
                                1.29e-03
                                </td>
                                <td style="text-align:left;">
                                -30%
                                </td>
                                <td style="text-align:left;">
                                0.886
                                </td>
                                <td style="text-align:left;">
                                -0.9%
                                </td>
                                <td style="text-align:left;">
                                3.82e-04
                                </td>
                                <td style="text-align:left;">
                                -73%
                                </td>
                                <td style="text-align:left;">
                                0.960
                                </td>
                                <td style="text-align:left;">
                                +0.3%
                                </td>
                                <td style="text-align:left;">
                                1.37e-03
                                </td>
                                <td style="text-align:left;">
                                -39%
                                </td>
                                </tr>
                                <tr>
                                <td style="text-align:left;">
                                MSR_P
                                </td>
                                <td style="text-align:right;">
                                50
                                </td>
                                <td style="text-align:left;">
                                0.938
                                </td>
                                <td style="text-align:left;">
                                +0.4%
                                </td>
                                <td style="text-align:left;">
                                4.15e-04
                                </td>
                                <td style="text-align:left;">
                                -78%
                                </td>
                                <td style="text-align:left;">
                                0.845
                                </td>
                                <td style="text-align:left;">
                                -6%
                                </td>
                                <td style="text-align:left;">
                                1.14e-04
                                </td>
                                <td style="text-align:left;">
                                -92%
                                </td>
                                <td style="text-align:left;">
                                0.957
                                </td>
                                <td style="text-align:left;">
                                +0.008%
                                </td>
                                <td style="text-align:left;">
                                8.11e-04
                                </td>
                                <td style="text-align:left;">
                                -64%
                                </td>
                                </tr>
                                <tr>
                                <td style="text-align:left;">
                                MSR_E
                                </td>
                                <td style="text-align:right;">
                                50
                                </td>
                                <td style="text-align:left;">
                                0.936
                                </td>
                                <td style="text-align:left;">
                                +0.08%
                                </td>
                                <td style="text-align:left;">
                                1.17e-04
                                </td>
                                <td style="text-align:left;">
                                -94%
                                </td>
                                <td style="text-align:left;">
                                0.820
                                </td>
                                <td style="text-align:left;">
                                -8%
                                </td>
                                <td style="text-align:left;">
                                8.93e-05
                                </td>
                                <td style="text-align:left;">
                                -94%
                                </td>
                                <td style="text-align:left;">
                                0.954
                                </td>
                                <td style="text-align:left;">
                                -0.3%
                                </td>
                                <td style="text-align:left;">
                                0.00e+00
                                </td>
                                <td style="text-align:left;">
                                -100%
                                </td>
                                </tr>
                                </tbody>
                                </table>
                                …

# Supplementary Material

## Additional functions for sup-mat

``` r
get_overlap <- function(rend, tstart, rstart, tend) {
  
  if (rend < tstart || rstart > tend) {
    return (0)
  }
  
  innerEnd <- min(rend, tend)
  outerEnd <- max(rend, tend)
  innerStart <- min(rstart, tstart)
  outerStart <- max(rstart, tstart)
  return (
    (innerEnd - innerStart) / (outerEnd - outerStart)
  )
}

parse_paf <- function(filepath, funcName, names) {
  gc()
  # raw <- filter_secondary_alignments(read_paf(filepath))
  raw <- filter_secondary_alignments(read_paf(filepath))
  
  df <- dplyr::as_tibble(raw) %>%
    mutate(rstart=as.integer(str_split(qname, "!", simplify=TRUE)[,3])) %>%
    mutate(rend=as.integer(str_split(qname, "!", simplify=TRUE)[,4])) %>%
    mutate(rtarget=str_split(qname, "!", simplify=TRUE)[,2]) %>%
    mutate(rstrand=str_split(qname, "!", simplify=TRUE)[,5]) %>% 
    mutate(rid=str_split(qname, "!", simplify=TRUE)[,1]) %>%
    rowwise() %>%
    mutate(overlap=get_overlap(rend, tstart, rstart, tend))
  
  df <- df %>% 
    mutate(correct=overlap>=0.1 && tname == rtarget) %>% 
    mutate(reduction=funcName) %>%
    mutate(chromosome=names[[rtarget]])
  gc()
  return(df)
  
}

parse_processed <- function(filepath, funcName, names) {
  df <- read_delim(filepath, delim="\t") %>%
    mutate(correct=as.logical(correct), reduction=funcName, chromosome=as.character(names[rtarget]))
  
  return(df)
}
```

Data loading

``` r
namesDic <- fromJSON(file="./data/case_study/names_chr.json")
mapping_raw <- parse_processed("./data/case_study/raw.processed.tab", "raw", namesDic)
mapping_hpc <- parse_processed("./data/case_study/hpc.processed.tab", "HPC", namesDic)
mapping_msr_p <- parse_processed("./data/case_study/msr_p.processed.tab", "MSR P", namesDic)
mapping_msr_e <- parse_processed("./data/case_study/msr_e.processed.tab", "MSR E", namesDic)
mapping_msr_f <- parse_processed("./data/case_study/msr_f.processed.tab", "MSR F", namesDic)
```

Centromere and stalk positions were obtained from the following file:
<http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.1/cytobands/cytoBandIdeo.bigBed>

``` r
centromeres <- tribble(
  ~chromosome, ~start, ~end,
  "1", 121796048, 126300487,
  "2", 92333543, 94673023,
  "3", 91738002, 96415026,
  "4", 49705154, 55199795,
  "5", 47039134, 49596625,
  "6", 58286706, 61058390,
  "7", 60414372, 63714499,
  "8", 44215832, 46325080,
  "9", 44951775, 47582595,
  "10", 39633793, 41664589,
  "11", 51035789, 54450838,
  "12", 34620838, 37202490,
  "13", 15547593, 17498291,
  "14", 10092112, 12708411,
  "15", 16678794, 17694466,
  "16", 35848286, 37829521,
  "17", 23892419, 27486939,
  "18", 15965699, 20933550,
  "19", 25817676, 29768171,
  "20", 26925852, 29099655,
  "21", 10962853, 11306205,
  "22", 12788180, 15711065,
  "X", 57820107, 60927025
)

stalks <- tribble(
  ~chromosome, ~start, ~end,
  "13", 5751447, 9368750,
  "14", 2077628, 2840421,
  "15", 2484618, 4728636,
  "21", 3084882, 5633495,
  "22", 4770731, 5743502 
)

order <- c(
  "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
  "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"
  )

reduction_order <- c("A", "MSR F2", "MSR P", "MSR E", "HPC", "raw")
```

## Table B.1

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="1">
</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="4">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

mapq=60

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="4">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

mapq\>=50

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="4">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

any mapq

</div>

</th>
</tr>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="1">
</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

fraction

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

error

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

fraction

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

error

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

fraction

</div>

</th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">

<div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">

error

</div>

</th>
</tr>
<tr>
<th style="text-align:left;">
MSR
</th>
<th style="text-align:left;">
f60
</th>
<th style="text-align:left;">
Pf60
</th>
<th style="text-align:left;">
e60
</th>
<th style="text-align:left;">
Pe60
</th>
<th style="text-align:left;">
f50
</th>
<th style="text-align:left;">
Pf50
</th>
<th style="text-align:left;">
e50
</th>
<th style="text-align:left;">
Pe50
</th>
<th style="text-align:left;">
f0
</th>
<th style="text-align:left;">
Pf0
</th>
<th style="text-align:left;">
e0
</th>
<th style="text-align:left;">
Pe0
</th>
</tr>
</thead>
<tbody>
<tr grouplength="5">
<td colspan="13" style="border-bottom: 1px solid;">
<strong>Drosophila - minimap</strong>
</td>
</tr>
<tr>
<td style="text-align:left;padding-left: 2em;" indentlevel="1">
HPC
</td>
<td style="text-align:left;">
0.957
</td>
<td style="text-align:left;">

-   0%
    </td>
    <td style="text-align:left;">
    2.27e-03
    </td>
    <td style="text-align:left;">

    -   0%
        </td>
        <td style="text-align:left;">
        0.963
        </td>
        <td style="text-align:left;">

        -   0%
            </td>
            <td style="text-align:left;">
            2.34e-03
            </td>
            <td style="text-align:left;">

            -   0%
                </td>
                <td style="text-align:left;">
                0.998
                </td>
                <td style="text-align:left;">

                -   0%
                    </td>
                    <td style="text-align:left;">
                    1.48e-02
                    </td>
                    <td style="text-align:left;">

                    -   0%
                        </td>
                        </tr>
                        <tr>
                        <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                        MSR_E
                        </td>
                        <td style="text-align:left;">
                        0.946
                        </td>
                        <td style="text-align:left;">
                        -1%
                        </td>
                        <td style="text-align:left;">
                        0.00e+00
                        </td>
                        <td style="text-align:left;">
                        -100%
                        </td>
                        <td style="text-align:left;">
                        0.954
                        </td>
                        <td style="text-align:left;">
                        -0.9%
                        </td>
                        <td style="text-align:left;">
                        0.00e+00
                        </td>
                        <td style="text-align:left;">
                        -100%
                        </td>
                        <td style="text-align:left;">
                        0.998
                        </td>
                        <td style="text-align:left;">
                        +0.05%
                        </td>
                        <td style="text-align:left;">
                        1.53e-02
                        </td>
                        <td style="text-align:left;">

                        -   3%
                            </td>
                            </tr>
                            <tr>
                            <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                            MSR_F
                            </td>
                            <td style="text-align:left;">
                            0.952
                            </td>
                            <td style="text-align:left;">
                            -0.5%
                            </td>
                            <td style="text-align:left;">
                            1.18e-03
                            </td>
                            <td style="text-align:left;">
                            -48%
                            </td>
                            <td style="text-align:left;">
                            0.960
                            </td>
                            <td style="text-align:left;">
                            -0.3%
                            </td>
                            <td style="text-align:left;">
                            1.37e-03
                            </td>
                            <td style="text-align:left;">
                            -41%
                            </td>
                            <td style="text-align:left;">
                            0.998
                            </td>
                            <td style="text-align:left;">
                            +0.008%
                            </td>
                            <td style="text-align:left;">
                            1.36e-02
                            </td>
                            <td style="text-align:left;">
                            -8%
                            </td>
                            </tr>
                            <tr>
                            <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                            MSR_P
                            </td>
                            <td style="text-align:left;">
                            0.950
                            </td>
                            <td style="text-align:left;">
                            -0.8%
                            </td>
                            <td style="text-align:left;">
                            4.90e-04
                            </td>
                            <td style="text-align:left;">
                            -78%
                            </td>
                            <td style="text-align:left;">
                            0.957
                            </td>
                            <td style="text-align:left;">
                            -0.6%
                            </td>
                            <td style="text-align:left;">
                            8.11e-04
                            </td>
                            <td style="text-align:left;">
                            -65%
                            </td>
                            <td style="text-align:left;">
                            0.998
                            </td>
                            <td style="text-align:left;">
                            -0.01%
                            </td>
                            <td style="text-align:left;">
                            1.39e-02
                            </td>
                            <td style="text-align:left;">
                            -6%
                            </td>
                            </tr>
                            <tr>
                            <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                            raw
                            </td>
                            <td style="text-align:left;">
                            0.958
                            </td>
                            <td style="text-align:left;">
                            +0.04%
                            </td>
                            <td style="text-align:left;">
                            2.27e-03
                            </td>
                            <td style="text-align:left;">
                            -0.04%
                            </td>
                            <td style="text-align:left;">
                            0.962
                            </td>
                            <td style="text-align:left;">
                            -0.08%
                            </td>
                            <td style="text-align:left;">
                            2.34e-03
                            </td>
                            <td style="text-align:left;">
                            +0.08%
                            </td>
                            <td style="text-align:left;">
                            0.997
                            </td>
                            <td style="text-align:left;">
                            -0.07%
                            </td>
                            <td style="text-align:left;">
                            1.17e-02
                            </td>
                            <td style="text-align:left;">
                            -21%
                            </td>
                            </tr>
                            <tr grouplength="5">
                            <td colspan="13" style="border-bottom: 1px solid;">
                            <strong>Drosophila - winnowmap</strong>
                            </td>
                            </tr>
                            <tr>
                            <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                            HPC
                            </td>
                            <td style="text-align:left;">
                            0.923
                            </td>
                            <td style="text-align:left;">

                            -   0%
                                </td>
                                <td style="text-align:left;">
                                1.51e-03
                                </td>
                                <td style="text-align:left;">

                                -   0%
                                    </td>
                                    <td style="text-align:left;">
                                    0.930
                                    </td>
                                    <td style="text-align:left;">

                                    -   0%
                                        </td>
                                        <td style="text-align:left;">
                                        1.59e-03
                                        </td>
                                        <td style="text-align:left;">

                                        -   0%
                                            </td>
                                            <td style="text-align:left;">
                                            0.989
                                            </td>
                                            <td style="text-align:left;">

                                            -   0%
                                                </td>
                                                <td style="text-align:left;">
                                                1.50e-02
                                                </td>
                                                <td style="text-align:left;">

                                                -   0%
                                                    </td>
                                                    </tr>
                                                    <tr>
                                                    <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                    MSR_E
                                                    </td>
                                                    <td style="text-align:left;">
                                                    0.905
                                                    </td>
                                                    <td style="text-align:left;">
                                                    -2%
                                                    </td>
                                                    <td style="text-align:left;">
                                                    1.42e-03
                                                    </td>
                                                    <td style="text-align:left;">
                                                    -6%
                                                    </td>
                                                    <td style="text-align:left;">
                                                    0.912
                                                    </td>
                                                    <td style="text-align:left;">
                                                    -2%
                                                    </td>
                                                    <td style="text-align:left;">
                                                    1.49e-03
                                                    </td>
                                                    <td style="text-align:left;">
                                                    -6%
                                                    </td>
                                                    <td style="text-align:left;">
                                                    0.983
                                                    </td>
                                                    <td style="text-align:left;">
                                                    -0.6%
                                                    </td>
                                                    <td style="text-align:left;">
                                                    1.44e-02
                                                    </td>
                                                    <td style="text-align:left;">
                                                    -4%
                                                    </td>
                                                    </tr>
                                                    <tr>
                                                    <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                    MSR_F
                                                    </td>
                                                    <td style="text-align:left;">
                                                    0.918
                                                    </td>
                                                    <td style="text-align:left;">
                                                    -0.6%
                                                    </td>
                                                    <td style="text-align:left;">
                                                    1.27e-03
                                                    </td>
                                                    <td style="text-align:left;">
                                                    -16%
                                                    </td>
                                                    <td style="text-align:left;">
                                                    0.925
                                                    </td>
                                                    <td style="text-align:left;">
                                                    -0.5%
                                                    </td>
                                                    <td style="text-align:left;">
                                                    1.30e-03
                                                    </td>
                                                    <td style="text-align:left;">
                                                    -18%
                                                    </td>
                                                    <td style="text-align:left;">
                                                    0.987
                                                    </td>
                                                    <td style="text-align:left;">
                                                    -0.2%
                                                    </td>
                                                    <td style="text-align:left;">
                                                    1.37e-02
                                                    </td>
                                                    <td style="text-align:left;">
                                                    -9%
                                                    </td>
                                                    </tr>
                                                    <tr>
                                                    <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                    MSR_P
                                                    </td>
                                                    <td style="text-align:left;">
                                                    0.905
                                                    </td>
                                                    <td style="text-align:left;">
                                                    -2%
                                                    </td>
                                                    <td style="text-align:left;">
                                                    1.33e-03
                                                    </td>
                                                    <td style="text-align:left;">
                                                    -12%
                                                    </td>
                                                    <td style="text-align:left;">
                                                    0.912
                                                    </td>
                                                    <td style="text-align:left;">
                                                    -2%
                                                    </td>
                                                    <td style="text-align:left;">
                                                    1.53e-03
                                                    </td>
                                                    <td style="text-align:left;">
                                                    -3%
                                                    </td>
                                                    <td style="text-align:left;">
                                                    0.983
                                                    </td>
                                                    <td style="text-align:left;">
                                                    -0.6%
                                                    </td>
                                                    <td style="text-align:left;">
                                                    1.40e-02
                                                    </td>
                                                    <td style="text-align:left;">
                                                    -7%
                                                    </td>
                                                    </tr>
                                                    <tr>
                                                    <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                    raw
                                                    </td>
                                                    <td style="text-align:left;">
                                                    0.949
                                                    </td>
                                                    <td style="text-align:left;">

                                                    -   3%
                                                        </td>
                                                        <td style="text-align:left;">
                                                        1.92e-03
                                                        </td>
                                                        <td style="text-align:left;">

                                                        -   27%
                                                            </td>
                                                            <td style="text-align:left;">
                                                            0.954
                                                            </td>
                                                            <td style="text-align:left;">

                                                            -   3%
                                                                </td>
                                                                <td style="text-align:left;">
                                                                1.99e-03
                                                                </td>
                                                                <td style="text-align:left;">

                                                                -   26%
                                                                    </td>
                                                                    <td style="text-align:left;">
                                                                    0.995
                                                                    </td>
                                                                    <td style="text-align:left;">
                                                                    +0.6%
                                                                    </td>
                                                                    <td style="text-align:left;">
                                                                    1.33e-02
                                                                    </td>
                                                                    <td style="text-align:left;">
                                                                    -12%
                                                                    </td>
                                                                    </tr>
                                                                    <tr grouplength="5">
                                                                    <td colspan="13" style="border-bottom: 1px solid;">
                                                                    <strong>Tandemtools -
                                                                    minimap</strong>
                                                                    </td>
                                                                    </tr>
                                                                    <tr>
                                                                    <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                    HPC
                                                                    </td>
                                                                    <td style="text-align:left;">
                                                                    0.870
                                                                    </td>
                                                                    <td style="text-align:left;">

                                                                    -   0%
                                                                        </td>
                                                                        <td style="text-align:left;">
                                                                        1.36e-03
                                                                        </td>
                                                                        <td style="text-align:left;">

                                                                        -   0%
                                                                            </td>
                                                                            <td style="text-align:left;">
                                                                            0.964
                                                                            </td>
                                                                            <td style="text-align:left;">

                                                                            -   0%
                                                                                </td>
                                                                                <td style="text-align:left;">
                                                                                1.56e-03
                                                                                </td>
                                                                                <td style="text-align:left;">

                                                                                -   0%
                                                                                    </td>
                                                                                    <td style="text-align:left;">
                                                                                    1.000
                                                                                    </td>
                                                                                    <td style="text-align:left;">

                                                                                    -   0%
                                                                                        </td>
                                                                                        <td style="text-align:left;">
                                                                                        9.00e-03
                                                                                        </td>
                                                                                        <td style="text-align:left;">

                                                                                        -   0%
                                                                                            </td>
                                                                                            </tr>
                                                                                            <tr>
                                                                                            <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                            MSR_E
                                                                                            </td>
                                                                                            <td style="text-align:left;">
                                                                                            0.885
                                                                                            </td>
                                                                                            <td style="text-align:left;">

                                                                                            -   2%
                                                                                                </td>
                                                                                                <td style="text-align:left;">
                                                                                                3.39e-03
                                                                                                </td>
                                                                                                <td style="text-align:left;">
                                                                                                +149%
                                                                                                </td>
                                                                                                <td style="text-align:left;">
                                                                                                0.962
                                                                                                </td>
                                                                                                <td style="text-align:left;">
                                                                                                -0.1%
                                                                                                </td>
                                                                                                <td style="text-align:left;">
                                                                                                3.53e-03
                                                                                                </td>
                                                                                                <td style="text-align:left;">
                                                                                                +127%
                                                                                                </td>
                                                                                                <td style="text-align:left;">
                                                                                                1.000
                                                                                                </td>
                                                                                                <td style="text-align:left;">

                                                                                                -   0%
                                                                                                    </td>
                                                                                                    <td style="text-align:left;">
                                                                                                    1.20e-02
                                                                                                    </td>
                                                                                                    <td style="text-align:left;">

                                                                                                    -   33%
                                                                                                        </td>
                                                                                                        </tr>
                                                                                                        <tr>
                                                                                                        <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                        MSR_F
                                                                                                        </td>
                                                                                                        <td style="text-align:left;">
                                                                                                        0.850
                                                                                                        </td>
                                                                                                        <td style="text-align:left;">
                                                                                                        -2%
                                                                                                        </td>
                                                                                                        <td style="text-align:left;">
                                                                                                        2.04e-03
                                                                                                        </td>
                                                                                                        <td style="text-align:left;">

                                                                                                        -   50%
                                                                                                            </td>
                                                                                                            <td style="text-align:left;">
                                                                                                            0.968
                                                                                                            </td>
                                                                                                            <td style="text-align:left;">
                                                                                                            +0.5%
                                                                                                            </td>
                                                                                                            <td style="text-align:left;">
                                                                                                            2.12e-03
                                                                                                            </td>
                                                                                                            <td style="text-align:left;">

                                                                                                            -   36%
                                                                                                                </td>
                                                                                                                <td style="text-align:left;">
                                                                                                                1.000
                                                                                                                </td>
                                                                                                                <td style="text-align:left;">

                                                                                                                -   0%
                                                                                                                    </td>
                                                                                                                    <td style="text-align:left;">
                                                                                                                    6.63e-03
                                                                                                                    </td>
                                                                                                                    <td style="text-align:left;">
                                                                                                                    -26%
                                                                                                                    </td>
                                                                                                                    </tr>
                                                                                                                    <tr>
                                                                                                                    <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                    MSR_P
                                                                                                                    </td>
                                                                                                                    <td style="text-align:left;">
                                                                                                                    0.898
                                                                                                                    </td>
                                                                                                                    <td style="text-align:left;">

                                                                                                                    -   3%
                                                                                                                        </td>
                                                                                                                        <td style="text-align:left;">
                                                                                                                        1.58e-03
                                                                                                                        </td>
                                                                                                                        <td style="text-align:left;">

                                                                                                                        -   16%
                                                                                                                            </td>
                                                                                                                            <td style="text-align:left;">
                                                                                                                            0.968
                                                                                                                            </td>
                                                                                                                            <td style="text-align:left;">
                                                                                                                            +0.4%
                                                                                                                            </td>
                                                                                                                            <td style="text-align:left;">
                                                                                                                            1.79e-03
                                                                                                                            </td>
                                                                                                                            <td style="text-align:left;">

                                                                                                                            -   15%
                                                                                                                                </td>
                                                                                                                                <td style="text-align:left;">
                                                                                                                                1.000
                                                                                                                                </td>
                                                                                                                                <td style="text-align:left;">

                                                                                                                                -   0%
                                                                                                                                    </td>
                                                                                                                                    <td style="text-align:left;">
                                                                                                                                    9.78e-03
                                                                                                                                    </td>
                                                                                                                                    <td style="text-align:left;">

                                                                                                                                    -   9%
                                                                                                                                        </td>
                                                                                                                                        </tr>
                                                                                                                                        <tr>
                                                                                                                                        <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                        raw
                                                                                                                                        </td>
                                                                                                                                        <td style="text-align:left;">
                                                                                                                                        0.936
                                                                                                                                        </td>
                                                                                                                                        <td style="text-align:left;">

                                                                                                                                        -   8%
                                                                                                                                            </td>
                                                                                                                                            <td style="text-align:left;">
                                                                                                                                            1.86e-03
                                                                                                                                            </td>
                                                                                                                                            <td style="text-align:left;">

                                                                                                                                            -   36%
                                                                                                                                                </td>
                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                0.984
                                                                                                                                                </td>
                                                                                                                                                <td style="text-align:left;">

                                                                                                                                                -   2%
                                                                                                                                                    </td>
                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                    2.09e-03
                                                                                                                                                    </td>
                                                                                                                                                    <td style="text-align:left;">

                                                                                                                                                    -   34%
                                                                                                                                                        </td>
                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                        1.000
                                                                                                                                                        </td>
                                                                                                                                                        <td style="text-align:left;">

                                                                                                                                                        -   0%
                                                                                                                                                            </td>
                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                            4.50e-03
                                                                                                                                                            </td>
                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                            -50%
                                                                                                                                                            </td>
                                                                                                                                                            </tr>
                                                                                                                                                            <tr grouplength="5">
                                                                                                                                                            <td colspan="13" style="border-bottom: 1px solid;">
                                                                                                                                                            <strong>Tandemtools -
                                                                                                                                                            winnowmap</strong>
                                                                                                                                                            </td>
                                                                                                                                                            </tr>
                                                                                                                                                            <tr>
                                                                                                                                                            <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                            HPC
                                                                                                                                                            </td>
                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                            0.775
                                                                                                                                                            </td>
                                                                                                                                                            <td style="text-align:left;">

                                                                                                                                                            -   0%
                                                                                                                                                                </td>
                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                1.32e-03
                                                                                                                                                                </td>
                                                                                                                                                                <td style="text-align:left;">

                                                                                                                                                                -   0%
                                                                                                                                                                    </td>
                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                    0.822
                                                                                                                                                                    </td>
                                                                                                                                                                    <td style="text-align:left;">

                                                                                                                                                                    -   0%
                                                                                                                                                                        </td>
                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                        1.82e-03
                                                                                                                                                                        </td>
                                                                                                                                                                        <td style="text-align:left;">

                                                                                                                                                                        -   0%
                                                                                                                                                                            </td>
                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                            0.997
                                                                                                                                                                            </td>
                                                                                                                                                                            <td style="text-align:left;">

                                                                                                                                                                            -   0%
                                                                                                                                                                                </td>
                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                8.37e-02
                                                                                                                                                                                </td>
                                                                                                                                                                                <td style="text-align:left;">

                                                                                                                                                                                -   0%
                                                                                                                                                                                    </td>
                                                                                                                                                                                    </tr>
                                                                                                                                                                                    <tr>
                                                                                                                                                                                    <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                    MSR_E
                                                                                                                                                                                    </td>
                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                    0.795
                                                                                                                                                                                    </td>
                                                                                                                                                                                    <td style="text-align:left;">

                                                                                                                                                                                    -   2%
                                                                                                                                                                                        </td>
                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                        2.28e-03
                                                                                                                                                                                        </td>
                                                                                                                                                                                        <td style="text-align:left;">

                                                                                                                                                                                        -   73%
                                                                                                                                                                                            </td>
                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                            0.846
                                                                                                                                                                                            </td>
                                                                                                                                                                                            <td style="text-align:left;">

                                                                                                                                                                                            -   3%
                                                                                                                                                                                                </td>
                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                2.52e-03
                                                                                                                                                                                                </td>
                                                                                                                                                                                                <td style="text-align:left;">

                                                                                                                                                                                                -   38%
                                                                                                                                                                                                    </td>
                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                    0.997
                                                                                                                                                                                                    </td>
                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                    -0.02%
                                                                                                                                                                                                    </td>
                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                    6.96e-02
                                                                                                                                                                                                    </td>
                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                    -17%
                                                                                                                                                                                                    </td>
                                                                                                                                                                                                    </tr>
                                                                                                                                                                                                    <tr>
                                                                                                                                                                                                    <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                                    MSR_F
                                                                                                                                                                                                    </td>
                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                    0.820
                                                                                                                                                                                                    </td>
                                                                                                                                                                                                    <td style="text-align:left;">

                                                                                                                                                                                                    -   6%
                                                                                                                                                                                                        </td>
                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                        1.83e-03
                                                                                                                                                                                                        </td>
                                                                                                                                                                                                        <td style="text-align:left;">

                                                                                                                                                                                                        -   38%
                                                                                                                                                                                                            </td>
                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                            0.867
                                                                                                                                                                                                            </td>
                                                                                                                                                                                                            <td style="text-align:left;">

                                                                                                                                                                                                            -   6%
                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                2.27e-03
                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                <td style="text-align:left;">

                                                                                                                                                                                                                -   25%
                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                    0.997
                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                    -0.06%
                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                    5.97e-02
                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                    -29%
                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                    </tr>
                                                                                                                                                                                                                    <tr>
                                                                                                                                                                                                                    <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                                                    MSR_P
                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                    0.780
                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                    +0.6%
                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                    1.62e-03
                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                    <td style="text-align:left;">

                                                                                                                                                                                                                    -   22%
                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                        0.829
                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                        +0.9%
                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                        2.09e-03
                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                        <td style="text-align:left;">

                                                                                                                                                                                                                        -   15%
                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                            0.997
                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                            -0.02%
                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                            8.65e-02
                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                            <td style="text-align:left;">

                                                                                                                                                                                                                            -   3%
                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                </tr>
                                                                                                                                                                                                                                <tr>
                                                                                                                                                                                                                                <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                                                                raw
                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                0.850
                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                <td style="text-align:left;">

                                                                                                                                                                                                                                -   10%
                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                    2.04e-03
                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                    <td style="text-align:left;">

                                                                                                                                                                                                                                    -   54%
                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                        0.890
                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                        <td style="text-align:left;">

                                                                                                                                                                                                                                        -   8%
                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                            1.95e-03
                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                            <td style="text-align:left;">

                                                                                                                                                                                                                                            -   7%
                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                0.999
                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                +0.2%
                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                4.60e-02
                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                -45%
                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                </tr>
                                                                                                                                                                                                                                                <tr grouplength="5">
                                                                                                                                                                                                                                                <td colspan="13" style="border-bottom: 1px solid;">
                                                                                                                                                                                                                                                <strong>Whole
                                                                                                                                                                                                                                                Human
                                                                                                                                                                                                                                                genome -
                                                                                                                                                                                                                                                minimap</strong>
                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                </tr>
                                                                                                                                                                                                                                                <tr>
                                                                                                                                                                                                                                                <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                                                                                HPC
                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                0.935
                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                <td style="text-align:left;">

                                                                                                                                                                                                                                                -   0%
                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                    1.85e-03
                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                    <td style="text-align:left;">

                                                                                                                                                                                                                                                    -   0%
                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                        0.942
                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                        <td style="text-align:left;">

                                                                                                                                                                                                                                                        -   0%
                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                            1.85e-03
                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                            <td style="text-align:left;">

                                                                                                                                                                                                                                                            -   0%
                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                1.000
                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                <td style="text-align:left;">

                                                                                                                                                                                                                                                                -   0%
                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                    1.46e-02
                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                    <td style="text-align:left;">

                                                                                                                                                                                                                                                                    -   0%
                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                        </tr>
                                                                                                                                                                                                                                                                        <tr>
                                                                                                                                                                                                                                                                        <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                                                                                                        MSR_E
                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                        0.926
                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                        -0.9%
                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                        6.92e-05
                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                        -96%
                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                        0.936
                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                        -0.7%
                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                        1.17e-04
                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                        -94%
                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                        0.999
                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                        -0.01%
                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                        1.76e-02
                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                        <td style="text-align:left;">

                                                                                                                                                                                                                                                                        -   20%
                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                            </tr>
                                                                                                                                                                                                                                                                            <tr>
                                                                                                                                                                                                                                                                            <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                                                                                                            MSR_F
                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                            0.930
                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                            -0.5%
                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                            1.09e-03
                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                            -41%
                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                            0.938
                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                            -0.4%
                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                            1.29e-03
                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                            -30%
                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                            1.000
                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                            -0.002%
                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                            1.51e-02
                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                            <td style="text-align:left;">

                                                                                                                                                                                                                                                                            -   4%
                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                </tr>
                                                                                                                                                                                                                                                                                <tr>
                                                                                                                                                                                                                                                                                <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                                                                                                                MSR_P
                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                0.929
                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                -0.6%
                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                2.20e-04
                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                -88%
                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                0.938
                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                -0.4%
                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                4.15e-04
                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                -78%
                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                0.999
                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                -0.02%
                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                1.55e-02
                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                <td style="text-align:left;">

                                                                                                                                                                                                                                                                                -   6%
                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                    </tr>
                                                                                                                                                                                                                                                                                    <tr>
                                                                                                                                                                                                                                                                                    <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                                                                                                                    raw
                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                    0.921
                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                    -1%
                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                    1.86e-03
                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                    +0.2%
                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                    0.927
                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                    -2%
                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                    1.86e-03
                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                    +0.6%
                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                    0.998
                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                    -0.2%
                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                    1.29e-02
                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                    -11%
                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                    </tr>
                                                                                                                                                                                                                                                                                    <tr grouplength="5">
                                                                                                                                                                                                                                                                                    <td colspan="13" style="border-bottom: 1px solid;">
                                                                                                                                                                                                                                                                                    <strong>Whole
                                                                                                                                                                                                                                                                                    Human
                                                                                                                                                                                                                                                                                    genome -
                                                                                                                                                                                                                                                                                    winnowmap</strong>
                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                    </tr>
                                                                                                                                                                                                                                                                                    <tr>
                                                                                                                                                                                                                                                                                    <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                                                                                                                    HPC
                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                    0.894
                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                    <td style="text-align:left;">

                                                                                                                                                                                                                                                                                    -   0%
                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                                        1.43e-03
                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                        <td style="text-align:left;">

                                                                                                                                                                                                                                                                                        -   0%
                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                            0.902
                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                            <td style="text-align:left;">

                                                                                                                                                                                                                                                                                            -   0%
                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                1.49e-03
                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                -   0%
                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                    0.988
                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                    -   0%
                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                        1.92e-02
                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                        -   0%
                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                            </tr>
                                                                                                                                                                                                                                                                                                            <tr>
                                                                                                                                                                                                                                                                                                            <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                                                                                                                                            MSR_E
                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                            0.795
                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                            -11%
                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                            6.33e-05
                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                            -96%
                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                            0.820
                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                            -9%
                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                            8.93e-05
                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                            -94%
                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                            0.971
                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                            -2%
                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                            2.08e-02
                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                            -   9%
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                </tr>
                                                                                                                                                                                                                                                                                                                <tr>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                                                                                                                                                MSR_F
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                0.874
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                -2%
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                2.81e-04
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                -80%
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                0.886
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                -2%
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                3.82e-04
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                -74%
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                0.984
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                -0.4%
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                1.94e-02
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                +0.9%
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                </tr>
                                                                                                                                                                                                                                                                                                                <tr>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                                                                                                                                                MSR_P
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                0.826
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                -8%
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                8.68e-05
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                -94%
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                0.845
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                -6%
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                1.14e-04
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                -92%
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                0.975
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                -1%
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                2.11e-02
                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                -   10%
                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                    </tr>
                                                                                                                                                                                                                                                                                                                    <tr>
                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                                                                                                                                                    raw
                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                    0.932
                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                    -   4%
                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                        1.75e-03
                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                        -   23%
                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                            0.937
                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                            -   4%
                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                1.79e-03
                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                -   20%
                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                    0.994
                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                    +0.7%
                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                    1.43e-02
                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                    -26%
                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                    </tr>
                                                                                                                                                                                                                                                                                                                                    <tr grouplength="5">
                                                                                                                                                                                                                                                                                                                                    <td colspan="13" style="border-bottom: 1px solid;">
                                                                                                                                                                                                                                                                                                                                    <strong>Whole
                                                                                                                                                                                                                                                                                                                                    Human
                                                                                                                                                                                                                                                                                                                                    genome
                                                                                                                                                                                                                                                                                                                                    (repeated
                                                                                                                                                                                                                                                                                                                                    regions) -
                                                                                                                                                                                                                                                                                                                                    minimap</strong>
                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                    </tr>
                                                                                                                                                                                                                                                                                                                                    <tr>
                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                                                                                                                                                                    HPC
                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                    0.619
                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                    -   0%
                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                        3.29e-04
                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                        -   0%
                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                            0.656
                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                            -   0%
                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                3.10e-04
                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                                -   0%
                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                    0.998
                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                                    -   0%
                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                        7.79e-02
                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                                        -   0%
                                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                                            </tr>
                                                                                                                                                                                                                                                                                                                                                            <tr>
                                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                                                                                                                                                                                            MSR_E
                                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                            0.618
                                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                            -0.2%
                                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                            1.41e-04
                                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                            -57%
                                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                            0.658
                                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                            +0.3%
                                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                            1.55e-04
                                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                            -50%
                                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                            0.997
                                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                            -0.03%
                                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                            8.23e-02
                                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                                            -   6%
                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                </tr>
                                                                                                                                                                                                                                                                                                                                                                <tr>
                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                                                                                                                                                                                                MSR_F
                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                0.601
                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                -3%
                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                2.18e-04
                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                -34%
                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                0.640
                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                -2%
                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                2.27e-04
                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                -27%
                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                0.998
                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                -0.02%
                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                8.15e-02
                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                                                -   5%
                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                    </tr>
                                                                                                                                                                                                                                                                                                                                                                    <tr>
                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                                                                                                                                                                                                    MSR_P
                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                    0.616
                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                    -0.5%
                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                    1.18e-04
                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                    -64%
                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                    0.656
                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                    +0.02%
                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                    1.99e-04
                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                    -36%
                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                    0.997
                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                    -0.08%
                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                    8.31e-02
                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                                                    -   7%
                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                        </tr>
                                                                                                                                                                                                                                                                                                                                                                        <tr>
                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                                                                                                                                                                                                        raw
                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                        0.514
                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                        -17%
                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                        1.98e-04
                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                        -40%
                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                        0.539
                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                        -18%
                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                        2.16e-04
                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                        -30%
                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                        0.981
                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                        -2%
                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                        6.69e-02
                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                        -14%
                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                        </tr>
                                                                                                                                                                                                                                                                                                                                                                        <tr grouplength="5">
                                                                                                                                                                                                                                                                                                                                                                        <td colspan="13" style="border-bottom: 1px solid;">
                                                                                                                                                                                                                                                                                                                                                                        <strong>Whole
                                                                                                                                                                                                                                                                                                                                                                        Human
                                                                                                                                                                                                                                                                                                                                                                        genome
                                                                                                                                                                                                                                                                                                                                                                        (repeated
                                                                                                                                                                                                                                                                                                                                                                        regions) -
                                                                                                                                                                                                                                                                                                                                                                        winnowmap</strong>
                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                        </tr>
                                                                                                                                                                                                                                                                                                                                                                        <tr>
                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                                                                                                                                                                                                        HPC
                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                        0.525
                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                                                        -   0%
                                                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                            1.24e-03
                                                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                                                            -   0%
                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                0.557
                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                                                                -   0%
                                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                    1.49e-03
                                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                                                                    -   0%
                                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                        0.950
                                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                                                                        -   0%
                                                                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                            1.19e-01
                                                                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                                                                            -   0%
                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                </tr>
                                                                                                                                                                                                                                                                                                                                                                                                <tr>
                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                                                                                                                                                                                                                                MSR_E
                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                0.366
                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                -30%
                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                6.35e-04
                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                -49%
                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                0.405
                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                -27%
                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                9.32e-04
                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                -37%
                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                0.911
                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                -4%
                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                1.38e-01
                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                                                                                -   17%
                                                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                                                    </tr>
                                                                                                                                                                                                                                                                                                                                                                                                    <tr>
                                                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                                                                                                                                                                                                                                    MSR_F
                                                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                    0.482
                                                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                    -8%
                                                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                    1.63e-03
                                                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                                                                                    -   31%
                                                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                        0.516
                                                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                        -7%
                                                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                        1.83e-03
                                                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                                                                                        -   23%
                                                                                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                            0.940
                                                                                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                            -1%
                                                                                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                            1.21e-01
                                                                                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                                                                                            -   2%
                                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                                </tr>
                                                                                                                                                                                                                                                                                                                                                                                                                <tr>
                                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                                                                                                                                                                                                                                                MSR_P
                                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                                0.415
                                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                                -21%
                                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                                9.45e-04
                                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                                -24%
                                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                                0.451
                                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                                -19%
                                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                                1.16e-03
                                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                                -22%
                                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                                0.920
                                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                                -3%
                                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                                1.39e-01
                                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                                                                                                -   17%
                                                                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                                                                    </tr>
                                                                                                                                                                                                                                                                                                                                                                                                                    <tr>
                                                                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;padding-left: 2em;" indentlevel="1">
                                                                                                                                                                                                                                                                                                                                                                                                                    raw
                                                                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                                    0.648
                                                                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                                                                                                    -   23%
                                                                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                                        1.26e-03
                                                                                                                                                                                                                                                                                                                                                                                                                        </td>
                                                                                                                                                                                                                                                                                                                                                                                                                        <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                                                                                                        -   1%
                                                                                                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                                            0.672
                                                                                                                                                                                                                                                                                                                                                                                                                            </td>
                                                                                                                                                                                                                                                                                                                                                                                                                            <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                                                                                                            -   21%
                                                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                                                1.49e-03
                                                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                                                +0.3%
                                                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                                                0.968
                                                                                                                                                                                                                                                                                                                                                                                                                                </td>
                                                                                                                                                                                                                                                                                                                                                                                                                                <td style="text-align:left;">

                                                                                                                                                                                                                                                                                                                                                                                                                                -   2%
                                                                                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                                                    8.09e-02
                                                                                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                                                                                    <td style="text-align:left;">
                                                                                                                                                                                                                                                                                                                                                                                                                                    -32%
                                                                                                                                                                                                                                                                                                                                                                                                                                    </td>
                                                                                                                                                                                                                                                                                                                                                                                                                                    </tr>
                                                                                                                                                                                                                                                                                                                                                                                                                                    </tbody>
                                                                                                                                                                                                                                                                                                                                                                                                                                    </table>

## Figure C.1

case study on high mapq reads

<img src="figures_files/figure-gfm/ridges_real-1.svg" width="100%" style="display: block; margin: auto;" />
The “A” entry is destined to be removed and is needed to scale
everything correctly between facets. To obtain the final version of the
figure the “A” track is removed using Inkscape. The remaining facets are
then moved closer to each other to avoid large gaps between facets. The
post-modification version of this plot can be seen
[here](./figures/figureC1.adapted.pdf).

## Figures D.1 to D.5

case study on all reads

### Figure D.1 (Raw)

<img src="figures_files/figure-gfm/sup_mat_raw-1.svg" width="100%" style="display: block; margin: auto;" />

### Figure D.2 (HPC)

<img src="figures_files/figure-gfm/sup_mat_hpc-1.svg" width="100%" style="display: block; margin: auto;" />

### Figure D.3 (MSR E)

<img src="figures_files/figure-gfm/sup_mat_msre-1.svg" width="100%" style="display: block; margin: auto;" />

### Figure D.4 (MSR P)

<img src="figures_files/figure-gfm/sup_mat_msrp-1.svg" width="100%" style="display: block; margin: auto;" />

### Figure D.5 (MSR F)

<img src="figures_files/figure-gfm/sup_mat_msrf-1.svg" width="100%" style="display: block; margin: auto;" />

## Figure E.1

<img src="figures_files/figure-gfm/plot_drosophila_ecoli_nanosim-1.png" style="display: block; margin: auto;" />
