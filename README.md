<h1 align="center">
  <br>
  <img src=https://user-images.githubusercontent.com/54723897/113879941-4e1af480-97bb-11eb-83f3-e0ec8772b7c4.gif width=500px>
  <br>    
</h1>

<h2 align="center">A benchmark dataset for cloud and cloud shadow semantic understanding</h2>

<p align="center">  
  • 
  <a href="#why-we-need-another-cloud-detection-dataset">Why we need another cloud detection dataset?</a> &nbsp;•  
  <a href="#characteristics">Characteristics</a> &nbsp;•
  <a href="#citation">Citation</a> &nbsp;•
  <a href="#credits">Credits</a>  
</p>

## Why we need another cloud detection dataset?

In the last years, we have witnessed the success of deep learning in semantic segmentation thanks to the availability of large-scale human-annotated datasets such as [**Cityscapes**](https://www.cityscapes-dataset.com/). However, in earth observation, more concretely in **cloud detection**, this revolution and remarkable gain in model efficiency have not happened yet.  Although many factors may be related to this, dataset availability is by far the central dilemma. To the best of our knowledge, current cloud detection datasets have some of these shortcomings:

- Regional spatial coverage.
- Lack of annotation richness.
- Insufficient metadata.
- Poor hand-crafted labeling.
- Lack of variety.
- Only designed to support standard supervised learning.
- Absence of a criteria of what we should to considerate as a cloud.

## Characteristics
	
- Provide 20 features, including SAR and multispectral data.
- Fifty thousand globally distributed image patches (511 x 511).
- Three different types of manual labeling: high-quality labels, scrabble labeling, and no labeling.
- Each image has the result of seven state-of-the-art cloud detection algorithm.
- Supports standard supervised learning and generalizations like weakly supervision, semi-supervised learning, few-shot learning, and novel approaches based on SAR-to-Optical image fusion.
- Freely available with love <3 using the STAC specification.
- 50 000 image patches (511 x 511) globally distributed.


## Metodology

<center>
  <br>
    <img src=https://user-images.githubusercontent.com/54723897/113933464-eb464f00-97f4-11eb-95a7-26ec235c47c9.png>
  <br>    
</center>

## input.npy description

<table><thead><tr><th></th><th>name</th><th>Index</th><th>Description</th><th>mean</th><th>sd</th></tr></thead><tbody><tr><td>1</td><td>B1</td><td>0</td><td></td><td>0.296554909974863</td><td>0.236695342071979</td></tr><tr><td>2</td><td>B2</td><td>1</td><td></td><td>0.282109446947271</td><td>0.245682208395963</td></tr><tr><td>3</td><td>B3</td><td>2</td><td></td><td>0.27402285332036</td><td>0.234278885807428</td></tr><tr><td>4</td><td>B4</td><td>3</td><td></td><td>0.297909889120112</td><td>0.249943283750605</td></tr><tr><td>5</td><td>B5</td><td>4</td><td></td><td>0.31032200778703</td><td>0.245354678243278</td></tr><tr><td>6</td><td>B6</td><td>5</td><td></td><td>0.342449812423938</td><td>0.23070012435287</td></tr><tr><td>7</td><td>B7</td><td>6</td><td></td><td>0.361181562353414</td><td>0.223915020757685</td></tr><tr><td>8</td><td>B8</td><td>7</td><td></td><td>0.346401739231663</td><td>0.216506874813553</td></tr><tr><td>9</td><td>B8A</td><td>8</td><td></td><td>0.370823464398451</td><td>0.216618439214885</td></tr><tr><td>10</td><td>B9</td><td>9</td><td></td><td>0.160060288466646</td><td>0.150134138981323</td></tr><tr><td>11</td><td>B10</td><td>10</td><td></td><td>0.0168750840165983</td><td>0.0332276294517062</td></tr><tr><td>12</td><td>B11</td><td>11</td><td></td><td>0.267365555783145</td><td>0.153045792042654</td></tr><tr><td>13</td><td>B12</td><td>12</td><td></td><td>0.20675507233452</td><td>0.129342459429327</td></tr><tr><td>14</td><td>VH</td><td>13</td><td></td><td>-22.878210105568</td><td>11.5174574548749</td></tr><tr><td>15</td><td>VV</td><td>14</td><td></td><td>-14.0824963873211</td><td>11.1807132211324</td></tr><tr><td>16</td><td>angle</td><td>15</td><td></td><td>38.5173328960839</td><td>4.77143493085221</td></tr><tr><td>17</td><td>CDI</td><td>16</td><td></td><td>-0.202116763846662</td><td>0.478734958132086</td></tr><tr><td>18</td><td>cloudshadow_direction</td><td>17</td><td></td><td>147.42553996465</td><td>48.5087154532189</td></tr><tr><td>19</td><td>elevation</td><td>18</td><td></td><td>1395.68532328893</td><td>1663.80277869266</td></tr><tr><td>20</td><td>Land Use</td><td>19</td><td></td><td>-</td><td>-</td></tr></tbody></table>

## manual.npy description

<table><thead><tr><th></th><th>Name</th><th>Description</th><th>Value</th></tr></thead><tbody><tr><td><br>1</td><td><br>Clear</td><td><br>All clear pixels, i.e. without cloud contamination or cloud shadows.</td><td><br>0</td></tr><tr><td><br>2</td><td><br>Thick Clouds</td><td><br>All cloudy pixels covered by thick clouds (does not include semi-transparent clouds or cloud shadows).</td><td><br>1</td></tr><tr><td><br>3</td><td><br>Thin Cloud</td><td><br>Clouds that are semi-transparent.</td><td><br>2</td></tr><tr><td><br>4</td><td><br>Cloud Shadows</td><td><br>All pixels contaminated by cloud shadows (not terrain shadows).</td><td><br>3</td></tr></tbody></table>


## Citation 

	COMMING SOON 
	
## Credits

	COMMING SOON 
