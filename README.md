# Ani_feces_roof
Repository for systematic review of animal fecal pathogens that have the potential to affect quality of roof harvested rain water.  Work conducted under ORISE and US EPA.

The Main Branch is where all the files that were directly held to creat the journal article.


N_Am_Meta_J-EPC0033278-QP-0_20230627_v03.csv
This is the R script that was used to conduct the metaanalysis.

Inputs
Campylobacter_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv
Coliform_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv
Cryptosporidium_STUDYN_DHD_J-EPC0033278-QP-1-0_220230627_v03.csv
EColi_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv
Enteriococci_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv
Giardia_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv
Salmonella_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv

These are the input files that were used used as part of the meta analysis.  The files are CSV tables that contain all the information that is needed to conduct the meta-analysis.  The files are seperated by pathogen.  The file naming convention is pathogen, STUDYN, initial of file creator (DHD PI of project), QAPP#, date, and version.

Outputs
Journal_forestplot.png
This is the forest plot that was created for the prevalence data for the joural article Fig. XX

Journal_EnumBoxplot.png
This is the boxplot that was created of the enumerated data for the journal article Fig. XX

Enumerated_box_display.csv
This is the table that is displayed in the main article and gives a numeric representation of the data in Journal_EnumBoxplot.png.

Enumerated_boxplot.csv
This table is the same as the information displayed in Enumerated_box_display.csv, but also includes data that could not fit in the journal article and includes the standard error and the names of the references for the outcomes that were part of the sysnesized data.

Total_Enumerated_Dataset.csv
This table is the same as Enumerated_boxplot.csv except for is also includes data for reports that were from the US and from groups that were from 1 study so were excluded from the table in the journal article but were included in the total.

Total_Prevalence_Dataset1.csv
This table contains all the information that was in Fig. XX Journal_forestplot.png, and also included the fixed effect model, random effect model confidence limits, and the references for all the reports that the data came from.
