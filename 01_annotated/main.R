


if (!exists("annotated_path", inherits = FALSE) || is.null(annotated_path)) {
  annotated_path <- "."
}


code_to_run <- list.files(annotated_path, pattern = "^\\d.*\\.R$")

for (code in code_to_run){
  source(code)
}

