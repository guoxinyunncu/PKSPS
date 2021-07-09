# PKSPS
PKSPS: a novel method for predicting kinase-specific phosphorylation sites combined with protein-protein interaction information and local sequence information.
# Requirement
```
numpy>=1.14.5
pandas>=1.0.3
networkx>=2.5
```
# Install 
python 3.6 in Linux and Windows.<br>
Since the program is written in python 3.6, python 3.6 with pip tool must be installed first. PKSPS uses the following dependencies: numpy, pandas, networkx, and random. You can install these packages first, by the following commands:<br>
```python
pip install numpy
pip install pandas
pip install networkx
pip install random
```
# Running PKSPS
Open cmd in windows or terminal in Linux, then download all the data and code in PKSPS-master to the local address, cd to where the code and data are stored and run:<br>
`
python predict.py
`
## Example:
```
python predict.py 
Please enter the substrate protein for inquiry:’TP53’
Please enter a sequence of queries:’ PSVEPPLsQETFSDL’
Please select a threshold for the output:0.01
```
Prediction results will show in the cmd or terminal

# Announcements
* Make sure the data and code are in one folder, or enter the exact data address when you run the code;<br>
* If an error occurs with the CMD runtime, consider running predict.py using the Python editor;<br>
* If you want to predict the catalytic kinase for the phosphorylation site, you need to provide the protein name of the substrate `Gene name`, and the sequence of amino acids around it `~-7~Site~+7~`;<br>
* The setting of the threshold can be selected by the user. We recommend that the user choose between the high threshold `0.01` and the low threshold `0.02`. User can also try the value between the two thresholds as needed;<br>
* The accepted amino acids are: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y, and a virtual amino acid O. If the protein fragments contain other amino acids, the program only will predict fragments which contain above-mentioned 21 amino acids.
