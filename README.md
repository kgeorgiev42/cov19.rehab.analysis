# Exploring Rehabilitation Pathways from Electronic Health Records of Patients with COVID-19
Supporting code for the workshop in _Process-Oriented Data Science for Healthcare (PODS4H)_, part of _ICPM (International Conference on Process Mining) 2023_.
The code shows how to develop sample cohorts of hospitalised patients with COVID-19 and generate summary features describing rehabilitation need from routine contacts with an AHP (Allied Healthcare Professional). The Process Mining functionality uses **bupaR** (fuzzy mining) to generate an event log from the timestamped contacts data and visualise patient interactions (during Wave 1 or Wave 2) using frequency-based **process maps**. Other functions allow additional analytics of the event log, such as plotting the precedence matrix, dotted chart, throughput and idle time over patient groups.
## Routine data
This is a case study that uses dummy data, based on **NHS Lothian** linked routine data from Electronic Health Records. In our case study, all data was collected from EHRs and national registries previously anonymised by the DataLoch service (**Edinburgh, United Kingdom**) and analysed in a Secure Data Environment. Details on these datasets are available on the [Metadata Catalogue](https://www.wiki.ed.ac.uk/display/DMCatalogue/2023.2%3A+DataLoch+Metadata+Catalogue+Navigation+Page). The provided dummy datasets describe:
  - The COVID-19 linked episodes of care (_episodes.rda_)
  - The performed lab tests for COVID-19 positivity (_lab_tests.rda_)
  - TrakCare pathways describing routine timestamped hospital visits (_trak_care.rda_)
  - AHP-related timestamped contact visits including session duration (_contacts.rda_)
  - TrakCare admission episodes (_trak_adm.rda_)
  - TrakCare discharge episodes (_trak_disch.rda_)
  - Scottish Morbidity Records (SMR-coded phenotypes), describing previous long-term conditions (_smr_phenotype.rda_)
  - Death records (_deaths.rda_)
## Requirements
This is an unofficial R library built in _RStudio and Rtools 4.3.0_ (recommended R >= 3.5.0). The required packages are listed in _DESCRIPTION_. Note that ```DiagrammeRsvg (>= 0.1), rsvg (>= 2.4.0)``` are typically required to produce pdf visualisations.
## Usage
<b>1.</b> Clone the project from this repository to your local directory using:
```
git clone https://github.com/JadeBlue96/Style-AI.git 
```
<b>2.</b> Install the R library using RStudio or a similar interpreter. This can be done by loading this repository as a project and clicking _Build -> Install and Restart_, or:
```
Rcmd.exe INSTALL --no-multiarch --with-keep.source cov19.rehab.analysis
```
<b>3.0.</b> The R library is now ready to use. The simplest way of launching the experiments is to use the pre-built _.rda_ files in the ```data/``` folder. However, if you would like to tweak the data, you can do so by modifying the parameters in ```generateData.R``` and running:
```
cov19.rehab.analysis::gen_cov_data()
```
to re-create all datasets or run the respective function of the modified dataset. Each R function contains adjustable parameters to change the target number of records, patients, records per patient, timestamps and likelihoods of sampling for categorical data.

<b>3.1.</b> To produce the cohort files (patient summary file and rehabilitation contacts file) you can run:
```
cov19.rehab.analysis::gen_cohorts()
```
This will re-run the code to sample COVID-19 positive patients and extract summary rehabilitation needs data based on _trak_care.rda_ and _contacts.rda_. This functionality can also be tweaked in ```generateEvalCohort.R```.

<b>3.2.</b> To produce the event log required to create the process maps you can run:
```
cov19.rehab.analysis::pm.generate_event_log(out_path='data/pm_event_log.rda', positive_only=1)
```
This will extract the timestamps from only the COVID-19 positive hospital episodes and save this in ```data/```.

<b>4.</b> The main functionalities include:
<b>4.1.</b> Generating a summary cohort table, stratified by COVID-19 wave.
```
cov19.rehab.analysis::pm.get_summary_table()
```
<b>4.2.</b> Plotting the frequent traces using a process map:
```
cov19.rehab.analysis::pm.plot_frequent_traces(out_path='results/freq_map.pdf', maf=5, minf=50, wave='Second')
```
The above will produce a process map describing traces within COVID-19 Wave 2, restricting the flows across the full map to at least 50 instances and requiring at least 5 patients with a specific interaction.

<b>4.3.</b> Plotting the pathways for a specific age group:
```
cov19.rehab.analysis::pm.plot_frequent_traces(out_path='results/freq_map_subage.pdf', maf=5, minf=50, wave='Second', age_group='80+')
```
The above will produce the same subset as the previous map but will show the traces present only in patients aged 80 and above.

<b>4.4.</b> Plot a precedence matrix for a specific wave:
```
cov19.rehab.analysis::pm.plot_precedence_matrix(out_path='results/prec_matrix.png', wave='Second')
```
The above will plot a table of antecedent-consequent pairs, representing the frequent interactions between specialist nodes, across those in Wave 2.

<b>4.5.</b> Other plots:
```
### Plot the total time required to complete the pathway across age groups
cov19.rehab.analysis::pm.plot_throughput_time(out_path='results/throughput_time.png', units='days') 
### Plot the cumulative time without an AHP-related session in the hospital
cov19.rehab.analysis::pm.plot_idle_time(out_path='results/idle_time.png', units='days')
### Plot a cumulative dotted chart showing the accumulated population over time by COVID-19 wave
cov19.rehab.analysis::pm.plot_dotted_chart(out_path='results/dotted_chart.pdf') 
```

## Event codes
AHP-related specialists:
- **NURSE** (care delivered by a specialist nurse)
- **PT** (physiotherapist)
- **OT** (occupational therapist)
- **DT** (dietician providing nutritional support)
- **SPL** (speech and language therapist)
- **CP** (clinical pharmacist)
- **PC** (specialist in palliative care)
- **INFSV** (infection prevention service)
-	**OTHER** (unidentified care provider)

## TODOs
- Provide documentation for R functions.
- Modify time to therapy and length of stay to provide more realistic pathways and admission episodes.
- Add functionality for linear regression analysis.
- Add coding and functionality for rehabilitation-related activities.

## Acknowledgements
We acknowledge the funding provided by the Sir Jules Thorn Charitable Trust PhD award (grant no. 21/01PhD).
