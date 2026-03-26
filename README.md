# Pulse_AnalyzeR

Example of applying pwa_plus.R using a for loop

```R
# load required functions from this repo
devtools::source_url("https://raw.githubusercontent.com/mkarmstrong/quick_R_plots/main/Splot.R")
devtools::source_url("https://raw.githubusercontent.com/mkarmstrong/quick_R_plots/main/Bplot.R")
devtools::source_url("https://raw.githubusercontent.com/mkarmstrong/quick_R_plots/main/BAplot.R")
devtools::source_url("https://raw.githubusercontent.com/mkarmstrong/quick_R_plots/main/BAplot.R")
devtools::source_url("https://raw.githubusercontent.com/mkarmstrong/quick_R_plots/main/BAplot.R")


# apply to multiple waveforms stored row wise in one .csv

# Load data
pwdata <- read.csv("C:/User/data/pulse_waves.csv")

# Create an empty data set to fill
pw_indices <- data.frame(matrix(data = NA, nrow = nrow(pwdata), ncol = 43))

# calculates the pw indcies and save them to pw_indices
for(i in 1:nrow(pwdata)) {
  # Extract ID
  id <- data.frame(ptid = pwdata[i, 1])
  print(unname(id))
  # Extract the waveform
  sig <- as.vector(na.omit(t(pwdata[i, -1])))
  # Calculate indices
  indices1 <- pwa_plus(sig, fs = 128, filt = T, norm = T)
  # Combined data
  results <- cbind(id, indices1)
  # Save into dataframe
  pw_indices[i, ] <- results[1,]
}

# col names
colnames(pw_indices) <- colnames(results)

# save data
write.csv(pw_indices, "R:/Data/KIDS/SCor Data/Katy_project/clean/pw_metrics.csv", row.names = F)
```
