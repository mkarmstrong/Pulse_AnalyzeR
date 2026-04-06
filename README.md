# Pulse_AnalyzeR

Example: apply the pulse wave analyzer using a for loop

```R
# install.packages("devtools") # if devtools not instaled, run once.

# load required functions from this repo
devtools::source_url("https://raw.githubusercontent.com/mkarmstrong/Pulse_AnalyzeR/refs/heads/main/pwa_functions.R")

# apply to multiple waveforms stored row wise in one .csv

# Load data (.csv) assumes ID variable is in first column with remaining columns containing the waveform
pwdata <- read.csv("R:/Data/KIDS/SCor Data/Katy_project/clean/pwaves_output.csv")

# Initialize data frame to fill
pw_indices <- data.frame(matrix(data = NA, nrow = nrow(pwdata), ncol = 44))

# Run analysis and save to pw_indices via a for loop (you could convert to lapply if prefered)
# Calculate pulse wave indices and save them to pw_indices via a for loop
for(i in 1:nrow(pwdata)) {
  id <- data.frame(ptid = pwdata[i, 1])       # Extract ID
  print(unname(id))                           # Print id to keep track
  sig <- as.vector(na.omit(t(pwdata[i, -1]))) # Extract the waveform
  indices1 <- pwa_plus(sig, 
                       ecgGated = F, 
                       fs = 200, 
                       filt = T, 
                       norm = T)              # Calculate indices (change defaults to suit, i.e. fs = 200/1000 etc.)
  results <- cbind(id, indices1)              # Combined data
  pw_indices[i, ] <- results[1,]              # Save into dataframe
}

# Get column names
colnames(pw_indices) <- colnames(results)

# Save data
write.csv(pw_indices, "R:/Data/KIDS/SCor Data/Katy_project/clean/pw_metrics.csv", row.names = F)
```
