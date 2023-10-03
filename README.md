# Ani_feces_roof Animal Taxonomy
Repository for systematic review of animal fecal pathogens that have the potential to affect quality of roof harvested rain water.  Work conducted under ORISE and US EPA.

These are all the files that were used to add the animal taxonomy to the literature review input tables in the main branch.
Campylobacter_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv  
Coliform_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv  
Cryptosporidium_STUDYN_DHD_J-EPC0033278-QP-1-0_220230627_v03.csv  
EColi_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv  
Enteriococci_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv  
Giardia_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv  
Salmonella_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv  


Animal_Names.R  
This is the R script that was used to add the animal taxonomy to the tables.  Using the script you can get most of the animal taxonomies to be added there are a few places where this doesn't work so make sure to got back through it once your done.

Animal Taxonomy DHD_J-EPC0033278-QP-1-0_20230407v03.csv
This table contains all the animal taxonomic informaiton that was used for this study.  If modifications are maded and more animal species are added then the they will have to be added to this table if the animal taxonomic information is going to be added using this script.

Input files
These files are pretty much the same as the input files for the main branch I just created seperate ones to not mess up the data in the main branch files while testing out the script.  There are only 4 input files for 4 pathogens because the other 3 were small data sets that were done by hand before it was decided to use the script to do it for me.

Campylobacter.csv  
Giardia.csv  
Salmonella.csv
Cryptosporidium.csv

Output files
These are a shorten version of the input files for this branch with a focus on adding the Animal taxonomic information.  The taxonomic information that was added to these output files was copied and pasted in to the main branch output files.

AddedAnimalNamesCampy.csv
AddedAnimalNamesCryp.csv
AddedAnimalNamesGia.csv
AddedAnimalNamesSal.csv
