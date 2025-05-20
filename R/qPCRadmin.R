############################################################################################################################
# Function to filter out all rows with the specified string(s) and write the result into a new file
# input is .txt or .csv results file from qPCR; only use Sample name(s) OR Target name(s), not both at once
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

    #check if .csv is from Taqman or SYBR run
    initial_rows <- readLines(input_file, n = 21)
    is_SYBR <- any(grepl('# Melt Stage Number:',initial_rows))

    if (is_SYBR) {
      qPCR_header <- read.csv(file = input_file,
                              header = TRUE, stringsAsFactors = FALSE, nrows = 20)
      qPCR_csv <- read.csv(file = input_file,
                           header = TRUE, stringsAsFactors = FALSE, skip = 21)
    } else {
      qPCR_header <- read.csv(file = input_file,
                              header = TRUE, stringsAsFactors = FALSE, nrows = 19)
      qPCR_csv <- read.csv(file = input_file,
                           header = TRUE, stringsAsFactors = FALSE, skip = 20)
    }
    # filter out all rows containing the specified strings in either the Sample or Target columns
    qPCR_csv <- qPCR_csv[!(qPCR_csv[['Sample']] %in% strings_to_remove | qPCR_csv[['Target']] %in% strings_to_remove),]
    write.csv(qPCR_header, output_file, row.names = FALSE)
    write.table(qPCR_csv, output_file,sep=",",append = TRUE, row.names = FALSE, col.names = TRUE)

  } else {
    # warning message if input file is neither .txt or .csv
    stop('ensure input_file is a .txt or .csv file')
  }
}

############################################################################################################################
# Function to filter rows with the specified string(s) and write the result into a new file
# input is .txt or .csv results file from qPCR; only use Sample name(s) OR Target name(s), not both at once
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

    #check if .csv is from Taqman or SYBR run
    initial_rows <- readLines(input_file, n = 21)
    is_SYBR <- any(grepl('# Melt Stage Number:',initial_rows))

    if (is_SYBR) {
      qPCR_header <- read.csv(file = input_file,
                              header = TRUE, stringsAsFactors = FALSE, nrows = 20)
      qPCR_csv <- read.csv(file = input_file,
                           header = TRUE, stringsAsFactors = FALSE, skip = 21)
    } else {
      qPCR_header <- read.csv(file = input_file,
                              header = TRUE, stringsAsFactors = FALSE, nrows = 19)
      qPCR_csv <- read.csv(file = input_file,
                           header = TRUE, stringsAsFactors = FALSE, skip = 20)
    }
    # filter to include only rows containing the specified strings in either the Sample or Target columns
    qPCR_csv <- qPCR_csv[qPCR_csv[['Sample']] %in% strings_to_include | qPCR_csv[['Target']] %in% strings_to_include,]
    write.csv(qPCR_header, output_file, row.names = FALSE)
    write.table(qPCR_csv, output_file,sep=",",append = TRUE, row.names = FALSE, col.names = TRUE)

  } else {
    # warning message if input file is neither .txt or .csv
    stop('ensure input_file is a .txt or .csv file')
  }
}

############################################################################################################################
# Function to change qPCR results.csv, .xls, or .xlsx file to a .txt file with only the minimum required data for quickPCR analysis
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
  cleaned_data <- data %>%
    filter(!grepl("NRT|NTC", .data[[sample_col]])) %>%
    select(Sample = all_of(sample_col), Target = all_of(target_col), Cq = all_of(cq_col))
  
  # Remove any additional information by dropping incomplete rows of the main data
  cleaned_data <- cleaned_data[complete.cases(cleaned_data), ]
  
  # Save the cleaned data as a .txt file
  write.table(cleaned_data, txt_output, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
}

############################################################################################################################
# Combines qPCR results.csv files into one .txt file
# files_to_combine and reporters_to_exclude should be input as vectors
combPCR <- function(files_to_combine, out_txt, reporters_to_exclude) {
  # set reporters_to_exclude as emtpy if not specified
  if(!hasArg(reporters_to_exclude)){reporters_to_exclude <- c()}

  #check if .csv is from Taqman or SYBR run
  initial_rows <- readLines(files_to_combine[1], n = 21)
  is_SYBR <- any(grepl('# Melt Stage Number:',initial_rows))

  # Initialize an empty list to store individual dataframes
  combined_data <- list()

  if (is_SYBR) {
    # Loop through each file
    for (file_path in files_to_combine) {
      # Check if the file exists before attempting to read it
      if (file.exists(file_path)) {
        # Read data
        qPCR_data <- read.csv(file_path, header = FALSE, stringsAsFactors = FALSE, skip = 21)
        # Extract the first non-skipped row as header
        headers <- qPCR_data[1, ]
        # Remove the first row from data as it is the header
        qPCR_data <- qPCR_data[-1, ]
        # Set the extracted headers to be the actual column names of the dataframe
        colnames(qPCR_data) <- headers
        # Define the column names to be removed
        columns_to_remove <- c("Well", "Well Position", "Omit", "Task", "Reporter", "Quencher", "Amp Status", "Amp Score",
                               "Curve Quality", "Result Quality Issues", "Cq Confidence", "Cq Mean", "Cq SD", "Auto Threshold",
                               "Threshold", "Auto Baseline", "Baseline Start", "Baseline End","Tm1","Tm2","Tm3","Tm4","Tm5")
        # Remove the specified columns by name
        qPCR_data <- qPCR_data[!(grepl("NRT", qPCR_data$'Sample') | grepl("NTC", qPCR_data$'Sample'))
                               , !(names(qPCR_data) %in% columns_to_remove)]
        # Remove rows with unwanted reporters
        if (length(reporters_to_exclude) > 0) {
          qPCR_data <- qPCR_data[!apply(qPCR_data, 1, function(row) any(row %in% reporters_to_exclude)), ]
        }
        # Append the edited data to the list
        combined_data[[length(combined_data) + 1]] <- qPCR_data
      } else {
        warning(paste("File not found:", file_path))
      }
    }
  } else {
    # Loop through each file
    for (file_path in files_to_combine) {
      # Check if the file exists before attempting to read it
      if (file.exists(file_path)) {
        # Read data
        qPCR_data <- read.csv(file_path, header = FALSE, stringsAsFactors = FALSE, skip = 20)
        # Extract the first non-skipped row as header
        headers <- qPCR_data[1, ]
        # Remove the first row from data as it is the header
        qPCR_data <- qPCR_data[-1, ]
        # Set the extracted headers to be the actual column names of the dataframe
        colnames(qPCR_data) <- headers
        # Define the column names to be removed
        columns_to_remove <- c("Well", "Well Position", "Omit", "Task", "Reporter", "Quencher", "Amp Status", "Amp Score",
                               "Curve Quality", "Result Quality Issues", "Cq Confidence", "Cq Mean", "Cq SD", "Auto Threshold",
                               "Threshold", "Auto Baseline", "Baseline Start", "Baseline End")
        # Remove the specified columns by name
        qPCR_data <- qPCR_data[!(grepl("NRT", qPCR_data$'Sample') | grepl("NTC", qPCR_data$'Sample'))
                               , !(names(qPCR_data) %in% columns_to_remove)]
        # Remove rows with unwanted reporters
        if (length(reporters_to_exclude) > 0) {
          qPCR_data <- qPCR_data[!apply(qPCR_data, 1, function(row) any(row %in% reporters_to_exclude)), ]
        }
        # Append the edited data to the list
        combined_data[[length(combined_data) + 1]] <- qPCR_data
      } else {
        warning(paste("File not found:", file_path))
      }
    }
  }
  # Combine all the edited data frames into one
  final_combined_data <- do.call(rbind, combined_data)
  # Write the combined data to a .txt file
  write.table(final_combined_data, out_txt, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}

#############################################################################################################################
# Adds a string (to use as an identifier) to each row of the 'Sample' column immediately proceeding sample name
# Reads in .txt file and outputs as new .txt file
IDtxt <- function(in_file, out_file, id) {
  # Read the combined .txt file into a dataframe
  combined_data <- read.table(in_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # Append character string to each value in the 'Sample' column
  combined_data$Sample <- paste0(combined_data$Sample, id)
  # Save the modified dataframe back to the .txt file
  write.table(combined_data, out_file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}

#############################################################################################################################
# Combine .txt files across experiments AFTER having added a unique character string in 'Sample' column
# If a unique ID is not needed, use combPCR to go straight from multiple .csv files to one .txt
combtxt <- function(txt_files_to_combine, out_file) {
  # Initialize an empty list to store dataframes
  txt_list <- list()
  # Read each .txt file and store the dataframe in the list
  for (file in txt_files_to_combine) {
    data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    txt_list <- append(txt_list, list(data))
  }
  # Combine all dataframes in the list into one dataframe
  combined_data <- do.call(rbind, txt_list)
  # Write the combined dataframe to a new .txt file
  write.table(combined_data, out_file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}
