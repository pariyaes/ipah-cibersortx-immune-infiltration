# Save package and R-version information for reproducibility.

if (!dir.exists("results")) dir.create("results", recursive = TRUE)

capture.output(
  sessionInfo(),
  file = file.path("results", "sessionInfo.txt")
)
