<h1 align="center">
  <br>
  <img src=https://user-images.githubusercontent.com/54723897/113879941-4e1af480-97bb-11eb-83f3-e0ec8772b7c4.gif width=500px>
  <br>    
</h1>

<h2 align="center">A global dataset for cloud and cloud shadow semantic understanding</h2>

<p align="center">  
  • 
  <a href="#why-we-need-another-cloud-detection-dataset">Introduction</a> &nbsp;•  
  <a href="#characteristics">Instructions</a> &nbsp;•
  <a href="#citation">Citation</a> &nbsp;•
  <a href="#credits">Credits</a>  
</p>

## Introduction

This repository provides the code necessary for selecting, downloading, and generating the cloudSEN12 image patches. The code is entirely reproducible; anybody interested in developing local or national versions of cloudSEN12 just has to make changes to the `data/cloudsen12_initial.csv` file. The [GH discussion section](https://github.com/cloudsen12/dataset/discussions) contains all our debates in the control quality stage.

## Instructions

1) Install all the libraries inside the `main.R` file.
2) Read the dataset (`data/cloudsen12_initial.csv`) generated by the CDE team after image tile selection.
3) Run `ip_creator` function. It will return all features of an image patch in cloudSEN12 as a list of stars objects: 
    - **s2l1c:** Sentinel-2 Level 1C data.
    - **s2l2a:** Sentinel-2 Level 2A data.
    - **extra:** Additional features described in Table 3 (See our paper [here]()).
    - **s1:** Sentinel-1.
6) Finally the metadata is generated by `metadata_creator` function.
7) Repeat 3 and 4 49400 times.

## Citation 

	COMMING SOON 
	
## Acknowledgment

This project gratefully acknowledges:

<img src=https://user-images.githubusercontent.com/16768318/153642319-9bb91ef6-a400-47ff-a080-9b4406390153.svg width=20%>

**for computing resources**

<img src=https://user-images.githubusercontent.com/16768318/153673173-e9069a03-daa7-4893-93ef-246248d48351.png width=20%>

**for rgee and rgeeExtra software**



