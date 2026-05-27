library(dplyr)
library(ggplot2)
library(data.table)
source('../../../lib.R')

# Simple case. Peak centered at position 50.
d1 <- data.table(seqnames = 'chrX', strand = '+',
                 start = c(46,  47,  48,  49,  50,  51,  52,  53,  54),
                 end   = c(90,  90,  90,  90,  90,  90,  90,  90,  90),
                 reads = c(1,    1,   5,   8,  10,   8,   5,   1,  1))

d2 <- standardize_positions(d1, side = "left", window = 8); d1$condition <- 'Unstandardized'; d2$condition <- 'Standardized'
d <- bind_rows(d2, d1); d$condition <- factor(d$condition, levels = c('Unstandardized', 'Standardized'))
p <- ggplot(d, aes(start, reads)) + 
     theme_minimal() + geom_col(color = 'black', fill = 'white') + 
     scale_x_continuous(breaks = seq(min(d1$start)-1, max(d1$start)+1, by = 1), limits = c(min(d1$start)-1, max(d1$start)+1)) + 
     facet_grid(~condition)
p



# Two closely spaced peaks centered at positions 48 and 53.
d1 <- data.table(seqnames = 'chrX', strand = '+',
                 start = c(46,  47,  48,  49,  50,    51,  52,  53,  54, 55),
                 end   = c(90,  90,  90,  90,  90,    90,  90,  90,  90, 90),
                 reads = c(1,    3,   10,  3,   1,     1,   3,  10,  3,   1))

d2 <- standardize_positions(d1, side = "left", window = 8); d1$condition <- 'Unstandardized'; d2$condition <- 'Standardized'
d <- bind_rows(d2, d1); d$condition <- factor(d$condition, levels = c('Unstandardized', 'Standardized'))
p <- ggplot(d, aes(start, reads)) + 
  theme_minimal() + geom_col(color = 'black', fill = 'white') + 
  scale_x_continuous(breaks = seq(min(d1$start)-1, max(d1$start)+1, by = 1), limits = c(min(d1$start)-1, max(d1$start)+1)) + 
  facet_grid(~condition)
p


# Three closely spaced peaks centered at positions 48, 53, and 58.
d1 <- data.table(seqnames = 'chrX', strand = '+',
                 start = c(46,  47,  48,  49,  50,    51,  52,  53,  54, 55,    56, 57, 58, 59, 60),
                 end   = c(90,  90,  90,  90,  90,    90,  90,  90,  90, 90,    90, 90, 90, 90, 90),
                 reads = c(1,    3,   10,  3,   1,     1,   3,  10,  3,   1,     1,  3,  10,  3,  1))

d2 <- standardize_positions(d1, side = "left", window = 8); d1$condition <- 'Unstandardized'; d2$condition <- 'Standardized'
d <- bind_rows(d2, d1); d$condition <- factor(d$condition, levels = c('Unstandardized', 'Standardized'))
p <- ggplot(d, aes(start, reads)) + 
  theme_minimal() + geom_col(color = 'black', fill = 'white') + 
  scale_x_continuous(breaks = seq(min(d1$start)-1, max(d1$start)+1, by = 1), limits = c(min(d1$start)-1, max(d1$start)+1)) + 
  facet_grid(~condition)
p


# Two peaks with one rising out of the tail of the taller peak at position 54.
d1 <- data.table(seqnames = 'chrX', strand = '+',
                 start = c(46,  47,  48,  49,  50,  51,  52,  53,  54,  55,     56),
                 end   = c(90,  90,  90,  90,  90,  90,  90,  90,  90,  90,    90),
                 reads = c(1,    3,   5,  10,  20,  10,   5,  3,    9,   4,     3))

d2 <- standardize_positions(d1, side = "left", window = 8); d1$condition <- 'Unstandardized'; d2$condition <- 'Standardized'
d <- bind_rows(d2, d1); d$condition <- factor(d$condition, levels = c('Unstandardized', 'Standardized'))
p <- ggplot(d, aes(start, reads)) + 
  theme_minimal() + geom_col(color = 'black', fill = 'white') + 
  scale_x_continuous(breaks = seq(min(d1$start)-1, max(d1$start)+1, by = 1), limits = c(min(d1$start)-1, max(d1$start)+1)) + 
  facet_grid(~condition)
p



# Chatter.
d1 <- data.table(seqnames = 'chrX', strand = '+',
                 start = c(46,  47,  48,  49,  50,  51,  52,  53,  54),
                 end   = c(90,  90,  90,  90,  90,  90,  90,  90,  90),
                 reads = c(1,    2,   1,  2,    1,   2,   1,  2,    1))

d2 <- standardize_positions(d1, side = "left", window = 8); d1$condition <- 'Unstandardized'; d2$condition <- 'Standardized'
d <- bind_rows(d2, d1); d$condition <- factor(d$condition, levels = c('Unstandardized', 'Standardized'))
p <- ggplot(d, aes(start, reads)) + 
  theme_minimal() + geom_col(color = 'black', fill = 'white') + 
  scale_x_continuous(breaks = seq(min(d1$start)-1, max(d1$start)+1, by = 1), limits = c(min(d1$start)-1, max(d1$start)+1)) + 
  facet_grid(~condition)
p



# Chatter followed by a tall peak at position 56.
d1 <- data.table(seqnames = 'chrX', strand = '+',
                 start = c(46,  47,  48,  49,  50,  51,  52,  53,  54, 55, 56, 57),
                 end   = c(90,  90,  90,  90,  90,  90,  90,  90,  90, 90, 90, 90),
                 reads = c(1,    2,   1,  2,    1,   2,   1,  2,    1,  5, 10,  5))

d2 <- standardize_positions(d1, side = "left", window = 8); d1$condition <- 'Unstandardized'; d2$condition <- 'Standardized'
d <- bind_rows(d2, d1); d$condition <- factor(d$condition, levels = c('Unstandardized', 'Standardized'))
p <- ggplot(d, aes(start, reads)) + 
  theme_minimal() + geom_col(color = 'black', fill = 'white') + 
  scale_x_continuous(breaks = seq(min(d1$start)-1, max(d1$start)+1, by = 1), limits = c(min(d1$start)-1, max(d1$start)+1)) + 
  facet_grid(~condition)
p




