


if (!exists("annotated_path", inherits = FALSE) || is.null(annotated_path)) {
  annotated_path <- "."
}


code_to_run <- list.files(annotated_path, pattern = "^\\d.*\\.R$", full.names = TRUE)

for (code in code_to_run){
  message("Running: ", code)
  source(code)
}

message("Finished main.R")

