# Program Name: test_processing_code
# Creation Date: 08/24/2021
# Full name of author: William (Jimmy) Fulp
# Project or Protocol: Caskey 820
# Purpose or description of program: Process B cell test data
# Location of program:https://github.com/FredHutch/Caskey820Analysis
# Location of input data: networks/cavd/studies/cvd856/test_dat


library(data.table)
library(tidyverse)
library(readxl)
library(CytoML)
library(ggcyto)
library(flowWorkspace)
library(flowCore)
library(assertthat)



# User Params -------------------------------------------------------------

# Should output results overwrite files
override <- FALSE

flow_path <- '/Volumes/networks/cavd/studies/cvd856/test_data'
# Need to exclude old manifest files
manifest_date <- '20211112-06'
output_path = '~/Temp/cvd856/flow_results'
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# Saving output and messages to log file
if (!interactive()) {
  # Adding timestamp to log so it stops overriding itself
  log_path <- file.path(output_path, paste0('log_', format(Sys.time(), "%Y%m%d_%H%M%S"),'.txt'))
  log_file <- file(log_path)
  sink(log_file,  split = T)
  sink(log_file,  type = 'message')
}


# Data Reading and Processing ---------------------------------------------

all_files <- list.files(path = flow_path, recursive = TRUE, full.names = TRUE) %>%
  str_subset('zip', negate = TRUE) %>%
  # remove temp files
  str_subset('~', negate = TRUE)

all_files_df <- tibble(
  path = all_files,
  file_name = basename(all_files)
) %>%
  mutate(
    upload_date = basename(dirname(all_files)),
    file_ext = file_name %>%
      str_sub(start = str_locate(file_name, '\\.[^\\.]*$')[,'start'] + 1)
  ) %>%
  dplyr::filter(
    file_name %>% str_detect('Controls|\\~|xml.xml', negate = TRUE),
    # dropping old manifests
    !(file_name %>% tolower() %>% str_detect('manifest') &
      upload_date != manifest_date)
  )



# manifest processing and checking ----------------------------------------
manifest_path <- all_files_df %>%
  dplyr::filter(file_name %>% tolower() %>% str_detect('manifest')) %>%
  pull(path)

manifest_donor <- read_excel(
  path =  manifest_path, sheet = 'Donor Files'
  ) %>%
  mutate(
    file_name = paste0(File, '.', FileExtension %>% tolower())
  )


manifest_xml <- read_excel(
  path = manifest_path, sheet = 'XML File'
) %>%
  mutate(
    file_name = paste0(File, '.', FileExtension %>% tolower())
  )


manifest <- bind_rows(manifest_donor, manifest_xml) %>%
  mutate(
    FileExtension = FileExtension %>% str_to_lower(),
    # creating tube and plate manually for now
    Tube = case_when(
      # for sort report if two 00X then first one is tube and second is plate
      FileType == 'Sort Report' & File %>% str_count('_00.') == 2 ~
        File %>% str_sub(start = str_locate(File, '_00._Sort')[,'start'] + 1,
                         end = str_locate(File, '_00._Sort')[,'start'] + 3),
      # for other files tube may be at end
      FileType != 'Sort Report' & File %>% str_detect('_00.$') ~
        File %>% str_sub(str_locate(File, '_00.$')[,'start'] + 1),
      FileType != 'Gating Coord.' ~ '001'
    ),
    Plate = case_when(
      # plate at end of sort report files
      FileType == 'Sort Report' ~
        File %>% str_sub(str_locate(File, '_00.$')[,'start'] + 1),
    )
  )



# linking manifest with files ---------------------------------------------


# printing extra files that are not in manifest
extra_manifest <- anti_join(manifest, all_files_df,
                            by = c('file_name', 'FileExtension' = 'file_ext'))
if (nrow(extra_manifest) > 0) {
  cat("The following files listed in the manifest but can't be found:\n")
  cat(paste0(extra_manifest$file_name, collapse = '\n'))
}

extra_files <- anti_join(all_files_df, manifest,
                         by = c('file_name', 'file_ext' = 'FileExtension')) %>%
  dplyr::filter(!file_name %>% tolower() %>% str_detect('manifest'))
if (nrow(extra_files) > 0) {
  cat('The following files are not listed in the manifest:\n')
  cat(paste0(extra_files$file_name, collapse = '\n'))
}

manifest_w_files <- inner_join(
  manifest,
  all_files_df, by = c('file_name', 'FileExtension' = 'file_ext')
  )




xml_files <- manifest_w_files %>%
  dplyr::filter(FileExtension == 'xml')
fcs_files <- manifest_w_files %>%
  dplyr::filter(FileExtension == 'fcs')
csv_counts_files <- manifest_w_files %>%
  dplyr::filter(FileExtension == 'csv',
                file_name %>% str_detect('Sort', negate = TRUE))
csv_sort_files <- manifest_w_files %>%
  dplyr::filter(FileExtension == 'csv',
                file_name %>% str_detect('Sort', negate = FALSE))




####### Processing by experiment--------------------------------------------------------------------
overall_start_time <- Sys.time()

for (i in xml_files$Experiment) {

  cat('\n\n\n')
  timestamp()
  start_time <- Sys.time()


  output_dir <- file.path(output_path, i)
  if (!override & dir.exists(output_dir) & length(dir(output_dir)) > 0) {
    cat('"overide" = F and output directory not empty, so skipping')
    next()
  }
  # might change later if exp name not enough info
  output_prefix <- i


  xml_file_here <- xml_files %>% filter(Experiment == i)
  fcs_files_here <- fcs_files %>% filter(Experiment == i)
  csv_counts_files_here <- csv_counts_files %>% filter(Experiment == i)
  csv_sort_files_here <- csv_sort_files %>% filter(Experiment == i)

  current_ptid <- unique(fcs_files_here$PTID)
  current_visit <- unique(fcs_files_here$Visit)

  cat(paste0("Processing Experiment: ", i,
             ",\n PTID(s): ", paste0(current_ptid, collapse = ', '),
             ',\n Visit: ', paste0(current_visit, collapse = ', '),
             '\n(', which(i == xml_files$Experiment),' of ',
             length(xml_files$Experiment), ')\n'))

  # initial checks
  if (n_distinct(dirname(fcs_files_here$path)) > 1) {
    cat('Error: All fcs files must be located in same directory folder\n')
    next
  }


  # Open diva worksheet -----------------------------------------------------
  # FlowJo-10.7.2

  cat('Opening Diva Workspace\n')

  diva_ws <- open_diva_xml(file = xml_file_here$path)

  fcs_names <- diva_get_sample_groups(diva_ws) %>%
    dplyr::filter(specimen %>% str_detect('G001 10X'))

  # testing all files in xml match files
  missing_fcs_files <- is.na(match(fcs_names$name, fcs_files_here$file_name))
  if (any(missing_fcs_files)) {
    cat('Warning: The following fcs files in xml file but not present:\n',
        paste0(fcs_names$name[missing_fcs_files], collapse = ', '),
        '\n', sep = '')
  }

  missing_xml_files <- is.na(match(fcs_files_here$file_name, fcs_names$name))
  if (any(missing_xml_files)) {
    cat('Warning: The following fcs files exist but not listed in xml experiment:\n',
        paste0(fcs_names$file_name[missing_xml_files], collapse = ', '),
        '\n', sep = '')
  }

  kw <-
    c(
      "$FIL",
      # "EXP_ASSAY_ID",
      # "Plate",
      # "PLATE ID",
      # "PLATE NAME",
      "TUBE NAME",
      "EXPERIMENT NAME",
      "Timepoint",
      "SampleID"
      # "Replicate",
      # "Comp",
      # "Type"
    )

    imported_workspace_global <-
    try(
      diva_to_gatingset(diva_ws,
                        name = 'G001 10X',
                        subset = fcs_files_here$file_name,
                        worksheet = 'global',
                        scale_level = 'gate',
                        verbose = TRUE,
                        truncate_max_range = FALSE)
  )
  if (inherits(imported_workspace_global, "try-error")) {
    cat('Error: Failed workspace import\n')
    next
  }

  redundant_nodes <- gs_check_redundant_nodes(imported_workspace_global) %>%
    map_dbl(length)

  if (sum(redundant_nodes) > 0) {
    cat('Error: Found redundant nodes in: ',
        paste0(names(redundant_nodes)[redundant_nodes > 0], collapse = ', '),
        '\n')
    next
  }

  # gs_pd <- pData(imported_workspace_global)
  #
  # channel_markers <- markernames(imported_workspace_global)
  #
  # gs_data <- gs_pop_get_data(obj = imported_workspace_global, y = "root")

  # gs_get_pop_paths(imported_workspace_global)


# Gating Plots ------------------------------------------------------------


  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Create plots for each fcs file to show the gating.
  cat("Creating gating plots:\n")
  plots_global <- list()
  for (s in sampleNames(imported_workspace_global)) {
    cat(s,'\n')
    suppressMessages(
      suppressWarnings(
        plots_global[[s]] <- autoplot(imported_workspace_global[[s]],
                                      bins = 64, strip.text = "gate")
      )
    )
    pdf(file = file.path(
      output_dir,
      paste0(s, "_globalWorksheet.pdf")
    ))
    suppressWarnings(
      print(plots_global[[s]])
    )
    dev.off()
  }


# Getting Counts ----------------------------------------------------------


  gs_counts <- gs_pop_get_count_fast(imported_workspace_global) %>%
    as_tibble() %>%
    mutate(Percent = Count / ParentCount * 100,
           Population = Population %>% str_trim(side = 'right')) %>%
    full_join(fcs_files_here,  by = c('name' = 'file_name'))


  cat('Writing out counts to ', output_dir, '\n')
  write_csv(gs_counts, file = file.path(output_dir, paste0(i, '_counts.csv')))




  # Read in Diva Counts ----------------------------------------------------------

  diva_csv_results <-
    map_dfr(csv_counts_files_here$path, function(x) {

      diva_counts_org <- as_tibble(fread(x, fill = TRUE, sep = ","))

      diva_counts <- diva_counts_org %>%
        filter(V3 != '', V1 != 'Population') %>%
        set_names(diva_counts_org %>% filter(V1 == 'Population')) %>%
        mutate(
          Diva_Count = as.numeric(`#Events`)
        ) %>%
        select(Diva_Population = Population,
               Diva_Parent = `Parent Name`,
               Diva_Count)

      left_join(
        diva_counts,
        diva_counts %>%
          select(Diva_Parent = Diva_Population,
                 Diva_ParentCount = Diva_Count),
        by = 'Diva_Parent'
      ) %>%
        mutate(
          Diva_Percent = Diva_Count / Diva_ParentCount * 100,
          path = x) %>%
        filter(Diva_Population != 'All Events')
    })
  diva_counts <- full_join(csv_counts_files_here, diva_csv_results, by = 'path')

  # Compare Counts to Lab ----------------------------------------------------------

  both_counts <- full_join(
    gs_counts %>%
      mutate(Population = basename(Population),
             Parent = ifelse(Parent %>% str_detect('root'),
                             'All Events', basename(Parent))) %>%
      select(Experiment, PTID, Visit, Tube, Population, Parent,
             Count, ParentCount, Percent),
    diva_counts %>%
      select(Experiment, PTID, Visit, Tube, Diva_Population, Diva_Parent,
             Diva_Count, Diva_ParentCount, Diva_Percent),
    by = c("Experiment", "PTID", "Visit", 'Tube',
           'Population' = 'Diva_Population', 'Parent' = 'Diva_Parent')
    ) %>%
    mutate(
      Count_Diff = abs(Count - Diva_Count),
      Count_Percent_Diff = Count_Diff / Count * 100
    )

  # Writing out mismatched parent/pops
  if (any(is.na(both_counts$Percent)))
    cat("The following populations are not available in the fcs files (Experiment_PTID_Visit_Tube): \n",
        both_counts %>%
          filter(is.na(Percent)) %>%
          transmute(paste0(Experiment, '_', PTID, '_', Visit, '_', Tube, ': ', Parent,'/',Population)) %>%
          distinct() %>% pull() %>% paste(collapse = '\n'), '\n'
    )
  if (any(is.na(both_counts$Diva_Percent)))
    cat("The following populations are not available in the diva csv files (Experiment_PTID_Visit_Tube): \n",
        both_counts %>%
          filter(is.na(Diva_Percent)) %>%
          transmute(paste0(Experiment, '_', PTID, '_', Visit, '_', Tube, ': ', Parent,'/',Population)) %>%
          distinct() %>% pull() %>% paste(collapse = '\n'), '\n'
    )



  # Getting spearman between diva and fcs
  spearman_label <-  both_counts %>%
    group_by(Experiment, PTID, Visit, Tube) %>%
    summarise(
      rho = cor(Percent, Diva_Percent),
      `.groups` = "drop"
    ) %>%
    mutate(rho_label = paste("rho = ", signif(rho,3)))

  both_counts_results <- full_join(
    both_counts,
    spearman_label ,
    by = c("Experiment", "PTID", "Visit", "Tube")
  )

  # Write out a CSV comparison with diva
  cat('Writing concordance with diva table to ', output_dir, '\n')
  write_csv(both_counts_results,
            file = file.path(
              output_dir,
              paste(output_prefix,
                    "concordance.csv",
                    sep = "_")
            ), na = '')



# Concordance plot with diva -----------------------------------------------------------------

  # Create a concordance plot of the extracted vs the DiVA proportions.
  cat('Generating VISC/DIVA concordance plot\n')

  pdf(file = file.path(
    output_dir,
    paste(output_prefix,
          "concordancePlot.pdf",
          sep = "_")
  ))

  concordance_plot <- both_counts_results %>%
    mutate(group = paste0(PTID, '_', Visit, '_', Tube, '(', rho_label, ')')) %>%
    filter(!is.na(Percent)) %>%
    ggplot(aes(x = Percent, y = Diva_Percent)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ylab("% of Parent Population (Diva)") +
    xlab("% of Parent Population (VISC)") +
    facet_wrap(~group) +
    # geom_text(data = spearman_label, aes(label = rho_label), y = 100, x = 15) +
    ggrepel::geom_label_repel(aes(label = Population), force = 5, size = 2.5, force_pull = 5)
    print(concordance_plot)

  dev.off()






  # Sort Report info --------------------------------------------------------


  sort_info <-
    map_dfr(csv_sort_files_here$path, function(x) {

      temp_sort <- as_tibble(fread(x, fill = TRUE, sep = ","))

      temp_out_vec <- c(
        'Threshold Count' = temp_sort$V2[temp_sort$V1 == 'Threshold Count'],
        'Processed Events Count(evt)' = temp_sort$V2[temp_sort$V1 == 'Processed Events Count(evt)'],
        'eOD-GT8++ Sorted' = temp_sort$V2[temp_sort$V2 %>% str_detect('eOD-GT8\\+\\+')] %>%
          str_replace('eOD-GT8\\+\\+ : ', ''),
        'KO11+ Sorted' = temp_sort$V3[temp_sort$V3 %>% str_detect('KO11\\+')] %>%
          str_replace('KO11\\+ : ', '')
      )

      temp_out_vec[temp_out_vec == 'Cont...'] <- 0
      temp_out <- as.numeric(temp_out_vec) %>% t() %>%
        as_tibble(.name_repair = 'minimal')
      names(temp_out) <- names(temp_out_vec)
      temp_out$path = x

      temp_out
    })

  sort_info_out <- full_join(csv_sort_files_here, sort_info, by = 'path')
  cat('Writing out sort reports to ', output_dir, '\n')
  write_csv(sort_info_out,
            file = file.path(output_dir, paste0(i, '_sort_reports.csv')),
            na = '')

}


cat('Overall Run Time:')
print(Sys.time() - overall_start_time)
