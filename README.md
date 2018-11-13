# Orthologs-Identification-Biopython

The ortholog scripts take as input, two different BLAST outputs in xml format and identify all pairs of putative orthologs from the two different species and generate the output in a text file. 

These scripts only consider e-values less than 10e-20.

The idea is simple, let’s say in Species A the proteins are noted A1, A2, A3 and in Species B they are B1, B2, B3. Let’s say A1 matches B1, B2, and B3 (i.e. "matches" means there is a HSP between them). Now let’s assume that B1 and B3 match A1. Pairs of putative orthologs would be:

A1 B1 

A1 B3

Note that A1 and B1 are not homologous pairs because A1 is not reciprocally a hit for B1. Now let’s assume that protein B1 from species B matches A3 from species A, then A3 is a homologous pair with B1. 

The output is a simple tab separated file with protein from one species in the first column and its ortholog pair from the other species in the second column. The protein definitions (e.g. 'gi| 671162305|ref|NP_007188796.2| CG30261, isoform I [Drosophila melanogaster]') are used to label the proteins. If a protein in species 1 does not have a pair in species 2, it is not printed out.

An optional final argument is included, such that if "FALSE" is passed as the last argument then the definition of a reciprocal match is no longer both a HSP from species A to B and species B to A, but is defined as a HSP from species A to B or species B to A.


## Prerequisites
Recommended Python 3.6.0 or later

Libraries used:

```
Biopython
```

## Installing
For Linux/BSD based systems, install the biopython package using pip.
```
pip install biopython
```
For windows systems, Anaconda or pip can be used.

The fasta files can be used to cross and generate the XML files using the `make_blastdb.sh` file.
Additional dependencies are required to use this shell script. 

## Running the tests

Run the script using the command line.
The `ortholog_dict.py` uses a dictionary data structure to efficiently find the putative ortholog pairs and <b>should be used</b>. The `ortholog_array.py` uses a simple array and set data structures to find the putative ortholog pairs which is not efficient and <b>should not be used</b>.

In Windows, open the Command Prompt then enter:

```
C:\Python36\python.exe C:\Users\...Path....\ortholog_dict.py cross_1.xml cross_2.xml orthologs_output.txt [FALSE]
```

In Unix, Linux and other BSD based systems, open bash shell then enter

```
chmod +x orthologs.py
python3 ortholog_dict.py cross_1.xml cross_2.xml orthologs_output.txt [FALSE]
```
The `[FALSE]` flag is optional for redefining the definition of homolog pairs.

## Built With

* [Python 3.6](https://www.python.org/downloads/release/python-360/) - The Programming tool used

## Versioning

Version tracked online with GitHub

## Authors

* **Samridha Shrestha**

## License

This project is licensed under the Apahce 2.0 License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Python open source libraries
* Biopython NCBI package and Biopython Developers
* NCBI https://www.ncbi.nlm.nih.gov/
