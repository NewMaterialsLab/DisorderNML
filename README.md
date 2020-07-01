# DisorderNML
DisorderNML is an application to explore the strcutural, electronic and optical properties of off-stoichiometric materials combined with SOD (https://github.com/gcmt-group/sod) and VASP outputs. First, you need to use SOD to create all possible atoms distribution in supercells. Then, VASP needs to be used to fully optimize the supercells (CalcXX/Opt) and calculate their optoelectronic properties (CalcXX/DOS).
 
usage: DisorderedLab.py [-h] [-m MODE] [-b RDF_BIN] [-x RDF_MAX]
                        [-A RDF_INDEXA] [-B RDF_INDEXB]
                        MasterFile Temp
positional arguments:
  MasterFile            File containing all the OUTSOD files.
  Temp                  Temperature. Units: Kelvin

optional arguments:
  -h, --help            show this help message and exit
  -m MODE, --mode MODE  Property to be calculated: dos, optics or rdf.
  -b RDF_BIN, --rdf_bin RDF_BIN
                        Bin width for the rdf. Units: Angstrom
  -x RDF_MAX, --rdf_max RDF_MAX
                        Maximun distance to be explored by the rdf. Units:
                        Angstrom
  -A RDF_INDEXA, --rdf_indexA RDF_INDEXA
                        Element 1 index for rdf
  -B RDF_INDEXB, --rdf_indexB RDF_INDEXB
                        Element 2 index for rdf

You can contact J. J. Plata for more information or reporting bugs at: jplata@us.es.
