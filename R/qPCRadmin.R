############################################################################################################################
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
