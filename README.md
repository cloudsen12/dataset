<h1 align="center">
  <br>
  <img src=https://user-images.githubusercontent.com/54723897/113879941-4e1af480-97bb-11eb-83f3-e0ec8772b7c4.gif width=500px>
  <br>    
</h1>

<h2 align="center">CloudSEN12 a global benchmark dataset for cloud semantic understanding</h2>

<p align="center">  
  • 
  <a href="#why-we-need-another-cloud-detection-dataset">Why we need another cloud detection dataset?</a> &nbsp;•  
  <a href="#metodology">Metodology</a> &nbsp;•
  <a href="#citation">Citation</a> &nbsp;•
  <a href="#credits">Credits</a>  
</p>

## Why we need another cloud detection dataset?

In the last years, we have witnessed the success of deep learning in semantic segmentation thanks to the availability of large-scale human-annotated datasets such as [**Cityscapes**](https://www.cityscapes-dataset.com/). However, in earth observation, more concretely in **cloud detection**, this revolution and remarkable gain in model efficiency have not happened yet.  Although many factors may be related to this, dataset availability is by far the central dilemma. In our best understanding, current cloud detection datasets experience some of these following deficiencies:

- Regional spatial coverage.
- Lack of annotation richness.
- Insufficient metadata.
- Poor hand-crafted labeling.
- Lack of variety.
- Only designed to support standard supervised learning.
- Absence of a criteria of what we should to considerate as a cloud.

## What is the main goal of cloudSEN12?


### Some characteristics

- 50 000 image patches (511 x 511) globally distributed.
- 

Exploiting the ready-to-use and freely data offering for Google Earth Engine we will create CloudSEN12, a new large dataset to advance deep learning research in cloud detection. CloudSEN12 will contain 25 000 globally distributed image patches (511x511 pixels) including (1) Sentinel-1 dual-pol synthetic aperture radar data, (2) Sentinel-2 multi-spectral data, (3) auxiliary remote sensed data, (4) hand-crafted labeling, and (5) the results from seven state-of-the-art cloud detection algorithms. The data set will be publicly available at https://console.cloud.google.com/storage/browser/cloudsen12 and contributions of novel methods are welcome to continuously providing a state-of-the-art benchmark dataset for cloud segmentation.



## Metodology

<center>
  <br>
    <img src=https://user-images.githubusercontent.com/54723897/113518868-36f4c080-9589-11eb-94f4-93ba4cab3c05.png>
  <br>    
</center>

## Citation 

## Credits

