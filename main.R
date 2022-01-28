# /*
# MIT License
#
# Copyright (c) [2022] [CloudSEN12 team]
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


# R packages --------------------------------------------------------------
library(reticulate)
library(lubridate)
library(dplyr)
library(readr)
library(rgee)
library(stars)
library(sf)
source("utils.R")

# Path with the cloudSEN12 dataset
CLOUDSEN12_PATH <- "/media/csaybar/Elements SE/cloudSEN12/"

# 1. Initialize Earth Engine ----------------------------------------------
ee_Initialize()


# 2. Load cloudsen12 initial dataset (after image tile selection) ---------
cloudsen12_init <- read.csv("data/cloudsen12_initial.csv") %>% 
  as_tibble()

# 3. Select an CloudSEN12 image patch (IP) --------------------------------
index <- 1
cloudsen12_ip <- cloudsen12_init[index,]

# 4. Create a cloudSEN12 IP. ----------------------------------------------
cloudsen12_ip_stars <- ip_creator(dataset = cloudsen12_ip)

# 5. Create metadata object for each IP. ----------------------------------
cloudsen12_ip_metadata <- metadata_creator(
  dataset = cloudsen12_ip,
  raster_ref = cloudsen12_ip_stars
)
