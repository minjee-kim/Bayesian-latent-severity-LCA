

# Handelman’s Dentistry Data (Handelman et al., 1986)
# ---------------------------------------------------
# This dataset is reconstructed from published frequency counts in:
# Handelman, S.L., Leverett, D.H., Espeland, M.A., & Curzon, J.A. 1986. 
# Clinical radiographic evaluation of sealed carious and sound tooth surfaces. 
# The Journal of the American Dental Association, 113(5), 751–754. 1516–1519.
#
# The data are publicly available in aggregated form in the original article.
# This script programmatically reconstructs an expanded binary dataset for
# research and reproducibility purposes only. Do not redistribute without citation.


patterns <- c("00000", "00001", "00010", "00011", "00100", "00101", "00110",
              "00111", "01000", "01001", "01010", "01011", "01100", "01101",
              "01110", "01111", "10000", "10001", "10010", "10011", "10100",
              "10101", "10110", "10111", "11000", "11001", "11010", "11011", 
              "11100", "11101", "11110", "11111")
freq     <- c(1880, 789, 43, 75, 23, 63, 8, 22, 188, 191, 17, 67, 15, 
              85, 8, 56, 22, 26, 6, 14, 1, 20, 2, 17, 2, 20, 6, 27, 3, 
              72, 1, 100)

expanded <- rep(patterns, freq)
k <- nchar(patterns[1])
data <- do.call(rbind, strsplit(expanded, ""))
data <- as.data.frame(apply(data, 2, as.numeric))