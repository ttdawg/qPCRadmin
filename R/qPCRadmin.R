######################################################################################################################################
# Function to test quality of qPCR run form its results.csv file.
# Checks for positives in negative controls NRT and NTC, their difference from samples values, and avg/med SD of technical replicates
# Outputs both a data frame within R and a .csv file
qPCRqa <- function(results_file_input, stats_output) {
  # Function to read data, regardless if it's csv or excel
  read_data <- function(file) {
    if (grepl("\\.csv$", file)) {
      data <- read.csv(file, header = FALSE, stringsAsFactors = FALSE)
    } else if (grepl("\\.xls$|\\.xlsx$", file)) {
      data <- read_excel(file, col_names = FALSE)
    } else {
      stop("Unsupported file format")
    }
    return(data)
  }
  
  # Read the file
  data <- read_data(results_file_input)
  
  # Identify the header row (where 'Well' is found in the first column)
  header_row <- which(data[, 1] == "Well")[1]
  
  # If the header row is not found, stop execution
  if (is.na(header_row)) {
    stop("Header row not found in the data.")
  }
  
  read_data_at_header <- function(file) {
    if (grepl("\\.csv$", file)) {
      data <- read.csv(file, header = FALSE, stringsAsFactors = FALSE, skip = header_row - 1)
    } else if (grepl("\\.xls$|\\.xlsx$", file)) {
      data <- read_excel(file, col_names = FALSE, skip = header_row - 1)
    } else {
      stop("Unsupported file format")
    }
    return(data)
  }
  
  data <- read_data_at_header(results_file_input)
  
  # Extract the first non-skipped row as header
  headers <- data[1, ]
  # Remove the first row from data as it is the header
  data <- data[-1, ]
  # Set the extracted headers to be the actual column names of the dataframe
  colnames(data) <- headers
  
  # Define possible column name sets
  sample_cols <- c("Sample", "Sample Name")
  target_cols <- c("Target", "Target Name")
  cq_cols <- c("Cq", "Cт")
  
  # Find actual column names present in the data
  sample_col <- sample_cols[sample_cols %in% colnames(data)]
  target_col <- target_cols[target_cols %in% colnames(data)]
  cq_col <- cq_cols[cq_cols %in% colnames(data)]
  
  # If any required columns aren't found, stop execution
  if (length(sample_col) == 0 || length(target_col) == 0 || length(cq_col) == 0) {
    stop("Required columns not found in the data.")
  }
  
  data <- data[, c(sample_col, target_col, cq_col)]
  colnames(data) <- c("Sample", "Target", "Cq")
  
  # Convert Cq values to numeric, treat non-numeric 'UNDETERMINED' as NA
  data_with_NAs <- data
  data_with_NAs$Cq <- as.numeric(as.character(data$Cq))
  
  # Remove any additional information by dropping incomplete rows of the main data
  data_NA_rmv <- data[complete.cases(data_with_NAs), ]
  
  # Identify targets where NRT/NTC have numbers
  all_nrt_ntc_data <- data_with_NAs[grepl("NRT|NTC", data_with_NAs$Sample) & !is.na(data$Cq), ]
  pos_nrt_ntc_data <- data_NA_rmv[grepl("NRT|NTC", data_NA_rmv$Sample) & !is.na(data$Cq), ]
  
  # Calculate metrics excluding NRT and NTC
  sample_data <- data_with_NAs[!grepl("NRT|NTC", data_with_NAs$Sample), ]
  
  # Calculate the mean and median of the standard deviations between technical replicates on the plate
  std_dev_df <- data_NA_rmv %>%
    group_by(Sample, Target) %>%
    summarize(std_dev = sd(Cq), .groups = 'drop')
  
  mean_std_dev <- mean(std_dev_df$std_dev, na.rm = TRUE)
  median_std_dev <- median(std_dev_df$std_dev, na.rm = TRUE)
  
  # Create the plate_sd_stats data frame
  plate_sd_stats <- data.frame(
    Target = 'Plate Stats',
    Plate_Avg_SD = mean_std_dev,
    Plate_Med_SD = median_std_dev
  )
  
  targets_summary_Cqs <- sample_data %>%
    group_by(Target) %>%
    summarise(
      avg_Cq = mean(Cq, na.rm = TRUE),
      med_Cq = median(Cq, na.rm = TRUE)
    )
  
  # Combine the 3 data frames
  unique_targets <- unique(pos_nrt_ntc_data$Target)
  filtered_summary_Cq <- targets_summary_Cqs[targets_summary_Cqs$Target %in% unique_targets, ]
  final_df <- data.frame(
    Target = unique_targets,
    NRT = NA,
    NTC = NA,  
    Plate_Avg_SD = NA,
    Plate_Med_SD = NA,
    stringsAsFactors = FALSE)
    
  for(i in 1:nrow(final_df)) {
    target <- final_df$Target[i]
    
    # Extract Cq values for NRT and NTC for the current target
    nrt_value <- pos_nrt_ntc_data$Cq[pos_nrt_ntc_data$Target == target & pos_nrt_ntc_data$Sample == "NRT"]
    ntc_value <- pos_nrt_ntc_data$Cq[pos_nrt_ntc_data$Target == target & pos_nrt_ntc_data$Sample == "NTC"]
    
    # Assign values to NRT and NTC columns
    final_df$NRT[i] <- ifelse(length(nrt_value) > 0, nrt_value, NA)
    final_df$NTC[i] <- ifelse(length(ntc_value) > 0, ntc_value, NA)
  }
  
  final_df <- merge(final_df, filtered_summary_Cq, by = "Target", all.x = TRUE)
  colnames(final_df)[colnames(final_df) == "avg_Cq"] <- "Sample_Avg"
  colnames(final_df)[colnames(final_df) == "med_Cq"] <- "Sample_Med"
  
  combined_df <- merge(final_df, plate_sd_stats, by = c("Target", "Plate_Avg_SD", "Plate_Med_SD"), all = TRUE)
  # Reorder columns
  combined_df <- combined_df[, c("Target", "NRT", "NTC", "Sample_Avg", "Sample_Med", "Plate_Avg_SD", "Plate_Med_SD")]
  
  # Reorder rows so that Plate Stats is last
  plate_stats_row <- combined_df[combined_df$Target == "Plate Stats", ]
  other_rows <- combined_df[combined_df$Target != "Plate Stats", ]
  
  final_ordered_df <- rbind(other_rows, plate_stats_row)
  
  write.csv(final_ordered_df, stats_output, row.names = FALSE)
  
  # Create a dataframe name by removing '.csv' from the `stats_output`
  df_name <- sub("\\.csv$", "", stats_output)
  
  # Assign the dataframe to the R environment with the desired name
  assign(df_name, final_ordered_df, envir = .GlobalEnv)
  
  return(final_ordered_df)
}

######################################################################################################################################
# Function to change qPCR results.csv to a .txt file with only the minimum required data for quickPCR analysis
minData <- function(results_file_input, txt_output) {
  # Function to read data, regardless if it's csv or excel
  read_data <- function(file) {
    if (grepl("\\.csv$", file)) {
      data <- read.csv(file, header = FALSE, stringsAsFactors = FALSE)
    } else if (grepl("\\.xls$|\\.xlsx$", file)) {
      data <- read_excel(file, col_names = FALSE)
    } else {
      stop("Unsupported file format")
    }
    return(data)
  }
  
  # Read the file
  data <- read_data(results_file_input)
  
  # Identify the header row (where 'Well' is found in the first column)
  header_row <- which(data[, 1] == "Well")[1]
  
  # If the header row is not found, stop execution
  if (is.na(header_row)) {
    stop("Header row not found in the data.")
  }
  
  read_data_at_header <- function(file) {
    if (grepl("\\.csv$", file)) {
      data <- read.csv(file, header = FALSE, stringsAsFactors = FALSE, skip = header_row - 1)
    } else if (grepl("\\.xls$|\\.xlsx$", file)) {
      data <- read_excel(file, col_names = FALSE, skip = header_row - 1)
    } else {
      stop("Unsupported file format")
    }
    return(data)
  }
  
  data <- read_data_at_header(results_file_input)
  
  # Extract the first non-skipped row as header
  headers <- data[1, ]
  # Remove the first row from data as it is the header
  data <- data[-1, ]
  # Set the extracted headers to be the actual column names of the dataframe
  colnames(data) <- headers
  
  # Define possible column name sets
  sample_cols <- c("Sample", "Sample Name")
  target_cols <- c("Target", "Target Name")
  cq_cols <- c("Cq", "Cт")
  
  # Find actual column names present in the data
  sample_col <- sample_cols[sample_cols %in% colnames(data)]
  target_col <- target_cols[target_cols %in% colnames(data)]
  cq_col <- cq_cols[cq_cols %in% colnames(data)]
  
  # If any required columns aren't found, stop execution
  if (length(sample_col) == 0 || length(target_col) == 0 || length(cq_col) == 0) {
    stop("Required columns not found in the data.")
  }
  
  # Filter out rows with unwanted samples and select unified columns
  cleaned_data <- data[!grepl("NRT|NTC", data[[sample_col]]), ]
  cleaned_data <- cleaned_data[, c(sample_col, target_col, cq_col)]
  colnames(cleaned_data) <- c("Sample", "Target", "Cq")
  
  # Remove any additional information by dropping incomplete rows of the main data
  cleaned_data <- cleaned_data[complete.cases(cleaned_data), ]
  
  # Save the cleaned data as a .txt or .csv file based on the extension of `txt_output`
  if (grepl("\\.csv$", txt_output)) {
    write.csv(cleaned_data, txt_output, row.names = FALSE)
  } else {
    write.table(cleaned_data, txt_output, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  }
  
}

############################################################################################################################
# Function to filter out all rows with the specified string(s) and write the result into a new file
# input is .txt or .csv results file from minData; only use Sample name(s) OR Target name(s), not both at once
rmvData <- function(input_file, output_file, strings_to_remove) {
  # Ensure strings_to_remove is always a character vector
  if (is.character(strings_to_remove)) {
    strings_to_remove <- as.character(strings_to_remove)
  } else {
    stop("strings_to_remove must be a character vector or string")
  }
  # Helper function to get file extension
  get_file_extension <- function(filename) {
    if (grepl("\\.", filename)) {
      return(tolower(sub("^.*\\.(.*)$", "\\1", filename)))
    }
    return("")
  }
  
  # Check if file is .txt or .csv
  file_extension <- get_file_extension(input_file)
  
  if (file_extension == "txt") {
    # Read the input file into a data frame
    data <- read.delim(input_file, stringsAsFactors = FALSE)
    # Remove rows containing any of the specified strings in any column
    cleaned_data <- data[!(data[['Sample']] %in% strings_to_remove | data[['Target']] %in% strings_to_remove), ]
    # Write the cleaned data to a new file
    write.table(cleaned_data, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
    
  } else if (file_extension == 'csv') {
    
    qPCR_csv <- read.csv(input_file, stringsAsFactors = FALSE)
    
    # filter out all rows containing the specified strings in either the Sample or Target columns
    qPCR_csv <- qPCR_csv[!(qPCR_csv[['Sample']] %in% strings_to_remove | qPCR_csv[['Target']] %in% strings_to_remove),]
    write.csv(qPCR_csv, output_file, row.names = FALSE)
    
  } else {
    # warning message if input file is neither .txt or .csv
    stop('ensure input_file is a .txt or .csv file')
  }
}


############################################################################################################################
# Function to filter rows with the specified string(s) and write the result into a new file
# input is .txt or .csv results file from minData; only use Sample name(s) OR Target name(s), not both at once
incldData <- function(input_file, output_file, strings_to_include) {
  # Ensure strings_to_remove is always a character vector
  if (is.character(strings_to_include)) {
    strings_to_include <- as.character(strings_to_include)
  } else {
    stop("strings_to_remove must be a character vector or string")
  }
  # Helper function to get file extension
  get_file_extension <- function(filename) {
    if (grepl("\\.", filename)) {
      return(tolower(sub("^.*\\.(.*)$", "\\1", filename)))
    }
    return("")
  }
  # Check if file is .txt or .csv
  file_extension <- get_file_extension(input_file)
  
  if (file_extension == "txt") {
    # Read the input file into a data frame
    data <- read.delim(input_file, stringsAsFactors = FALSE)
    # Remove rows containing any of the specified strings in any column
    cleaned_data <- data[data[['Sample']] %in% strings_to_include | data[['Target']] %in% strings_to_include, ]
    # Write the cleaned data to a new file
    write.table(cleaned_data, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  } else if (file_extension == 'csv') {
    
    qPCR_csv <- read.csv(input_file, stringsAsFactors = FALSE)
    
    # filter to include only rows containing the specified strings in either the Sample or Target columns
    qPCR_csv <- qPCR_csv[qPCR_csv[['Sample']] %in% strings_to_include | qPCR_csv[['Target']] %in% strings_to_include,]
    write.csv(qPCR_csv, output_file, row.names = FALSE)
    
  } else {
    # warning message if input file is neither .txt or .csv
    stop('ensure input_file is a .txt or .csv file')
  }
}

############################################################################################################################
# Combines minData files into one file.
# Input can be either .csv, .txt, or both. Output can be either .csv or .txt.
# files_to_combine should be input as a vector
combPCR <- function(files_to_combine, out_file) {
  
  # Initialize an empty list to store individual dataframes
  combined_data <- list()
  
  # Helper function to get file extension
  get_file_extension <- function(filename) {
    if (grepl("\\.", filename)) {
      return(tolower(sub("^.*\\.(.*)$", "\\1", filename)))
    }
    return("")
  }
  
  # Loop through each file
  for (file_path in files_to_combine) {
    
    file_extension <- get_file_extension(file_path)
    
    if (file_extension == 'csv') {
      # Read data
      qPCR_data <- read.csv(file_path, stringsAsFactors = FALSE)
      # Append the edited data to the list
      combined_data[[length(combined_data) + 1]] <- qPCR_data
    } else if (file_extension == 'txt') {
      qPCR_data <- read.delim(file_path, stringsAsFactors = FALSE)
      # Append the edited data to the list
      combined_data[[length(combined_data) + 1]] <- qPCR_data
    } else {
      stop("Unsupported file type as input: ", file_extension)
    }
  }
  
  # Combine all the edited data frames into one
  final_combined_data <- do.call(rbind, combined_data)
  
  # Save the combined data based on the extension of `out_file`
  output_extension <- get_file_extension(out_file)
  
  # Save the combined data as a .txt or .csv file based on the extension of `out_file`
  if (output_extension == "csv") {
    write.csv(final_combined_data, out_file, row.names = FALSE, quote = FALSE)
  } else if (output_extension == "txt") {
    write.table(final_combined_data, out_file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  } else {
    stop("Unsupported file type as output: ", output_extension)
  }
}

#############################################################################################################################
# Adds a string (to use as an identifier) to each row of the 'Sample' column immediately proceeding sample name
# Reads in .csv or .txt file and outputs as new .csv or .txt file
dataID <- function(in_file, out_file, id) {
  
  # Helper function to get file extension
  get_file_extension <- function(filename) {
    if (grepl("\\.", filename)) {
      return(tolower(sub("^.*\\.(.*)$", "\\1", filename)))
    }
    return("")
  }
  
  in_file_extension <- get_file_extension(in_file)
  out_file_extension <- get_file_extension(out_file)
  
  if (in_file_extension == 'csv') {
    data_to_id <- read.csv(in_file, stringsAsFactors = FALSE)
  } else if (in_file_extension == 'txt') {
    # Read the combined .txt file into a dataframe
    data_to_id <- read.delim(in_file, stringsAsFactors = FALSE)
  } else {
    stop("Unsupported file format (input): ", in_file_extension)
  }
  
  # Append character string to each value in the 'Sample' column
  data_to_id$Sample <- paste0(data_to_id$Sample, id)
  
  # Write modified dataframe to the output file based on out_file extension type
  if (out_file_extension == 'csv') {
    write.csv(data_to_id, out_file, row.names = FALSE, quote = FALSE)
  } else if (out_file_extension == 'txt') {
    write.table(data_to_id, out_file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  } else {
    stop("Unsupported file format (output): ", out_file_extension)
  }
}
