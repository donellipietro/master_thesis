cat.script_title <- function(title) {
  len <- nchar(title)
  spacer <- strrep("%", len)
  
  cat(paste("\n%%%", spacer, "%%%", sep = ""))
  cat(paste("\n%% ", title,  " %%", sep = ""))
  cat(paste("\n%%%", spacer, "%%%\n\n", sep = ""))
}

cat.section_title <- function(title) {
  len <- nchar(title)
  spacer <- strrep("|", len)
  
  cat(paste("\n|||", spacer, "|||", sep = ""))
  cat(paste("\n|| ", title,  " ||", sep = ""))
  cat(paste("\n|||", spacer, "|||\n\n", sep = ""))
}

save(cat.script_title, cat.section_title, file = "utils/functions/cat_utilities.RData")