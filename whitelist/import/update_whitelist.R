update_whitelist <- function(whitelist, combined_filter) {
  # Update the overall whitelist filter columns to include the putative variant filter information
  whitelist$Manual_Review <- (whitelist$Manual_Review & combined_filter != "FAIL") | (whitelist$Whitelist & combined_filter == "MANUAL_REVIEW")
  whitelist$Whitelist <- whitelist$Whitelist & combined_filter == "PASS"
  whitelist$Putative_Manual_Review <- (whitelist$Putative_Manual_Review & combined_filter != "FAIL") | (whitelist$Putative_Whitelist & combined_filter == "MANUAL_REVIEW")
  whitelist$Putative_Whitelist <- whitelist$Putative_Whitelist & combined_filter == "PASS"
  return(whitelist)
}