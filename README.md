The code is structured in different modules that are orchestrated by the main.py file and it is built to
work on multi-chromosomal vcf files. The remaining part of the code is split across the files listed below:

- Contains functions to parse and transform the original vcf content to return info about the input
file and tables in csv format
- Contains functions to adapt the csv files for plotting and generate png outputs
- Contains side functions
- Optional function to return a mutational signature report
The main.py file requests arguments from the user and parses them into the functions that are imported
from the other modules. The whole program can be run with default arguments using the command
below:

```bash
python main.py ../inputs/chr20.variants_20072022.vcf
```
