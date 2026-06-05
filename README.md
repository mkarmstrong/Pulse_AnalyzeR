# Pulse_AnalyzeR

Example: apply the pulse wave analyzer using a for loop

```R
# install.packages("devtools") # if devtools not instaled, run once.

# load required functions from this repo
devtools::source_url("https://raw.githubusercontent.com/mkarmstrong/Pulse_AnalyzeR/refs/heads/main/pwa_functions.R")

# apply to multiple waveforms stored row wise in one .csv (with the id in the first column)

# Load data (.csv) assumes ID variable is in first column with remaining columns containing the waveform
pwdata <- read.csv("R:/Data/path/pwaves.csv")

# Initialize list to fill
pw_indices <- vector("list", nrow(pwdata))

# Calculate pulse wave indices and save them to pw_indices via a for loop
for(i in 1:nrow(pwdata)) {
  id <- data.frame(ptid = pwdata[i, 1])                             # Extract ID
  cat("Processing", i, "of", nrow(pwdata), ":", pwdata[i, 1], "\n") # Print id to keep track
  sig <- as.vector(na.omit(t(pwdata[i, -1])))                       # Extract the waveform
  indices1 <- pwa_plus(sig,                                         # the pulse wave
                       ecgGated = FALSE,                            # Ensembled using ecg (set T), using P foot (set F)
                       fs = 200,                                    # Sample rate (Hz)
                       filt = T,                                    # Apply a low pass filter (T/F)
                       norm = F)                                    # Normalize to 0-1 range (T/F)
  pw_indices[[i]] <- cbind(id, indices1)                            # Save into pw_indices list
}

# Combine
pw_indices <- do.call(rbind, pw_indices)

# Save data
write.csv(pw_indices, "R:/Data/path/pwave_indices.csv", row.names = F)
```

Alternatively, you could apply the analysis to one waveform. This is often a useful sanity check.
```R
dat <- read.csv("R:/Data/path/pwaves.csv")
pw <- as.vector(na.omit(t(dat[10, -1]))); plot(pw) # here you can select the individual waves from pwaves.csv. Currently selecting the 10th wave.
pwa_plus(pw, ecgGated = F, filt = F, verbose = T) # here you can chage the defaults
```

Load pre-trained LightGBM models
```R
# load model trained with 0-1 calibration
model <- lgb.load("https://raw.githubusercontent.com/mkarmstrong/Pulse_AnalyzeR/refs/heads/main/cfpwv_lgbm_model.txt.R")
# predict cfPWV on new data set
predictions <- predict(model, as.matrix(new_data))

# load model trained with mean-diastolic calibration
model <- lgb.load("https://raw.githubusercontent.com/mkarmstrong/Pulse_AnalyzeR/refs/heads/main/cfpwv_lgbm_model_nn.txt.R")
# predict cfPWV on new data set
predictions <- predict(model, as.matrix(new_data))
```
