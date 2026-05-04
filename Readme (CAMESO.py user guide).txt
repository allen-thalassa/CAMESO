The cation norm (catanorm) and mesonorm are two normative mineralogy calculation schemes for igneous petrology; and they are crucial for specific classification diagrams and granitoid studies. 

The CAMESO.py (CAtion and MEsonorm SOftware), is a user-friendly, efficient Python-based program for calculating catanorms and mesonorms; and requires the Numpy, Pandas, and tqdm packages.


To use CAMESO.py:

i) Paste data into the comma-separated values (CSV) format file 'input_template.csv' and save;

ii) Open 'config.txt' and set 'method = 1' for catanorm calculation or 'method = 2' for mesonorm calculation. When the mesonorm procedure is chosen, the user has the option to compute norms with amphiboles by setting "path = 1"; otherwise, set "path = 2". Finally, save 'config.txt';

iii) Run 'CAMESO.py' under python environment; after several tens of seconds, a CSV format 'output.csv' file is created in the same directory;

iv) Open 'output.csv', edit, and copy the results into your own spreadsheet files.


