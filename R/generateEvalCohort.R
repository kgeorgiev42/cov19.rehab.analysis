gen_cohorts <- function() {
  ### Load dummy datasets
  set.seed(42)
  print('Loading dummy data..')
  load(file='data/trak_adm.rda')
  load(file='data/trak_disch.rda')
  load(file='data/smr_phenotype.rda')
  load(file='data/trak_care.rda')
  load(file='data/contacts.rda')
  load(file='data/episodes.rda')
  load(file='data/lab_tests.rda')
  load(file='data/deaths.rda')
  ### Cohort generation pipeline
  lab_full <- coh.get_index_admissions(admissions, discharges)
  lab_full <- coh.get_deaths(lab_full)
  lab_full <- coh.get_comorbidities(lab_full)
  coh.parse_contacts(lab_full)
  coh.adjust_distribution()

}

coh.get_index_admissions <- function(admissions, discharges, na_pos='Y') {
  #### Merge with TRAK
  print('Merging admissions and discharges..')
  trak_adm_ids <- unique(admissions$pid)
  lab_data <- lab_tests %>%
    select(c(pid, dataSource, collectionDate, testPositive))

  trak_adm_sel <- admissions %>%
    select(c(pid, episode, AdmissionDate, AdmissionTime)) %>%
    mutate(AdmissionDateTime = ymd_hm(paste(AdmissionDate, AdmissionTime)))

  trak_disch_sel <- discharges %>%
    select(c(pid, episode, DischargeDate, DischargeTime)) %>%
    mutate(DischargeDateTime = ymd_hm(paste(DischargeDate, DischargeTime)))

  trak_full <- left_join(trak_adm_sel, trak_disch_sel, by=c('pid', 'episode'))
  trak_full <- trak_full %>%
    select(c(pid, AdmissionDate, AdmissionDateTime, DischargeDate, DischargeDateTime, episode)) %>%
    na.omit() %>%
    #filter(AdmissionStatus %in% c('D', 'A')) %>%
    mutate(AdmissionDate = as.Date(AdmissionDate)) %>%
    mutate(DischargeDate = as.Date(DischargeDate))

  ### Handle overlapping admissions
  print('Merging repeating same-day admissions..')
  trak_sel <- trak_full %>%
    arrange(pid, AdmissionDate) %>%
    group_by(pid) %>%
    mutate(indx = c(0, cumsum(as.numeric(lead(AdmissionDate)) >
                                cummax(as.numeric(DischargeDate)))[-n()])) %>%
    group_by(pid,indx) %>%
    summarise(AdmissionDate=min(AdmissionDate), DischargeDate=max(DischargeDate),
              AdmissionDateTime=min(AdmissionDateTime), DischargeDateTime=max(DischargeDateTime), episode=first(episode)) %>%
    select(-indx)

  print(sprintf('Merged repeating admissions: %.0f out of %.0f', dim(trak_sel[0]), dim(trak_full[0]))[1])

  ### Merge with covid tests
  print('Getting COVID-19 positive admissions within 10 days of admission..')
  lab_data <- left_join(lab_data, trak_sel, by=c('pid'))
  lab_data <- lab_data[!is.na(lab_data$collectionDate),]
  #### Filter positive tests
  lab_data$testPositive[is.na(lab_data$testPositive)] <- na_pos
  ### Build COVID-19 positivity intervals
  admissionInterval10 <- interval(as.Date(lab_data$AdmissionDate) - days(10), as.Date(lab_data$DischargeDate))
  lab_data$collectionDate <- ymd(lab_data$collectionDate)

  lab_data <- lab_data %>%
    mutate(testPositive_10d = ifelse((collectionDate >= (as.Date(AdmissionDate) - days(10))) & (collectionDate < as.Date(DischargeDate)), 1, 0))

  ### Filter 10d lab tests
  lab_sel <- lab_data %>%
    filter(!((dataSource == 'Test by Lighthouse Laboratory') & (testPositive_10d == 0))) %>%
    select(c(pid, dataSource, collectionDate, testPositive_10d, episode, AdmissionDate, AdmissionDateTime, DischargeDate, DischargeDateTime)) %>%
    arrange(pid, AdmissionDateTime, desc(testPositive_10d)) %>%
    group_by(pid) %>%
    slice(1) %>%
    ungroup()

  #### Merge with baseline set of pts
  trak_pt_ids <- unique(lab_sel$pid)
  ep_sel <- episodes %>%
    filter((age_at_first_test >= 18) & (pid %in% trak_pt_ids))

  lab_full <- left_join(ep_sel, lab_sel, by='pid') %>%
    filter(!is.na(episode))

  ### Build COVID-19 wave flags
  lab_full <- lab_full %>%
    mutate(isFirstWave = ifelse((AdmissionDate >= as.Date("2020-03-01")) & (AdmissionDate <= as.Date("2020-05-31")), 1, 0)) %>%
    mutate(isSecondWave = ifelse((AdmissionDate >= as.Date("2020-09-01")) & (AdmissionDate <= as.Date("2021-03-31")), 1, 0))

  print(sprintf('Index episodes, extracted: %.0f out of %.0f', dim(lab_full[0]), dim(trak_sel[0]))[1])
  #save(lab_full, file='data/val.rda')
  lab_full
}

coh.get_deaths <- function(lab_full) {
  print('Merging death events and getting patients surviving hospital stay..')
  lab_full <- left_join(lab_full, deaths, by='pid')
  ### Died during stay flag
  lab_full <- lab_full %>%
    mutate(died_during_stay = ifelse(as.Date(DischargeDate, '%Y-%m-%d') >=
                                       as.Date(date_of_death, '%Y-%m-%d'), 1, 0)) %>%
    mutate(died_during_stay = ifelse(is.na(died_during_stay), 0, died_during_stay)) %>%
    filter(died_during_stay == 0)

  ### 1-year mortality
  lab_full$AdmissionDateTime <- lab_full$AdmissionDateTime[!is.na(lab_full$AdmissionDateTime)]
  lab_full$DischargeDateTime <- lab_full$DischargeDateTime[!is.na(lab_full$DischargeDateTime)]
  lab_full$AdmissionDate <- as.Date(lab_full$AdmissionDate[!is.na(lab_full$AdmissionDate)])
  lab_full$DischargeDate <- as.Date(lab_full$DischargeDate[!is.na(lab_full$DischargeDate)])
  mortality_interval <- interval(lab_full$AdmissionDate, lab_full$AdmissionDate + days(365))
  lab_full$date_of_death <- as.Date(ymd(lab_full$date_of_death))

  lab_full <- lab_full %>%
    mutate(mortality_1yr=case_when(date_of_death %within% mortality_interval ~1,
                                   TRUE ~0))

  print(sprintf('Patients extracted: %.0f', dim(lab_full[0]))[1])
  #save(lab_full, file='data/val.rda')
  lab_full
}

coh.get_comorbidities <- function(lab_full) {
  print('Extracting previous long-term conditions..')
  smr_phenotypes$EventDate <- as.Date(smr_phenotypes$EventDate[!is.na(smr_phenotypes$EventDate)])
  smr_phenotypes <- left_join(smr_phenotypes, lab_full[, c('pid', 'AdmissionDate')], by='pid')

  smr_sel <- smr_phenotypes %>%
    filter(!is.na(AdmissionDate)) %>%
    filter(EventDate < AdmissionDate) %>%
    select(c(pid, PhenotypeName))

  smr_sel$PhenotypeName <- gsub(",", "", as.character(smr_sel$PhenotypeName))

  smr_cms <- within(smr_sel, {PhenotypeName <- as.character(PhenotypeName);
  Comorbidities <- ave(PhenotypeName, pid, FUN=toString)})

  smr_cms <- smr_cms %>%
    select(c(pid, Comorbidities)) %>%
    group_by(pid) %>%
    slice(1)

  lab_full <- left_join(lab_full, smr_cms, by='pid')

  lab_full <- lab_full %>%
    mutate(n_morbid=str_split(Comorbidities, ',') %>% lengths) %>%
    mutate(n_morbid=ifelse(is.na(Comorbidities), 0, n_morbid)) %>%
    mutate(is_multimorbid=ifelse(n_morbid > 1, 1, 0)) %>%
    mutate(is_complex_multimorbid=ifelse(n_morbid > 3, 1, 0)) %>%
    mutate(cm_chd = ifelse(grepl('Coronary heart disease not otherwise specified', Comorbidities, fixed=TRUE), 1, 0)) %>%
    mutate(cm_stable_angina = ifelse(grepl('Stable angina', Comorbidities, fixed=TRUE), 1, 0)) %>%
    mutate(cm_unstable_angina = ifelse(grepl('Unstable Angina', Comorbidities, fixed=TRUE), 1, 0)) %>%
    mutate(cm_heart_failure = ifelse(grepl('Heart failure', Comorbidities, fixed=TRUE), 1, 0)) %>%
    mutate(cm_ic_haem = ifelse(grepl('Intracerebral haemorrhage', Comorbidities, fixed=TRUE), 1, 0)) %>%
    mutate(cm_ischaemic_stroke = ifelse(grepl('Ischaemic stroke', Comorbidities, fixed=TRUE), 1, 0)) %>%
    mutate(cm_stroke_nos = ifelse(grepl('Stroke NOS', Comorbidities, fixed=TRUE), 1, 0)) %>%
    mutate(cm_mi = ifelse(grepl('Myocardial infarction', Comorbidities, fixed=TRUE), 1, 0)) %>%
    mutate(cm_diabetes = ifelse(grepl('Diabetes', Comorbidities, fixed=TRUE), 1, 0)) %>%
    mutate(cm_obesity = ifelse(grepl('Obesity', Comorbidities, fixed=TRUE), 1, 0)) %>%
    mutate(cm_esrd = ifelse(grepl('End stage renal disease', Comorbidities, fixed=TRUE), 1, 0)) %>%
    mutate(cm_urinary_incontinence = ifelse(grepl('Urinary Incontinence', Comorbidities, fixed=TRUE), 1, 0)) %>%
    mutate(cm_hip_fracture = ifelse(grepl('Fracture of hip', Comorbidities, fixed=TRUE), 1, 0)) %>%
    mutate(cm_osteoporosis = ifelse(grepl('Osteoporosis', Comorbidities, fixed=TRUE), 1, 0)) %>%
    mutate(cm_parkinsons = ifelse(grepl('Parkinson\'s disease', Comorbidities, fixed=TRUE), 1, 0)) %>%
    mutate(cm_dementia = ifelse(grepl('Dementia', Comorbidities, fixed=TRUE), 1, 0)) %>%
    mutate(cm_delirium = ifelse(grepl('Delirium not induced by alcohol and other psychoactive substances', Comorbidities, fixed=TRUE), 1, 0)) %>%
    mutate(cm_depression = ifelse(grepl('Depression', Comorbidities, fixed=TRUE), 1, 0)) %>%
    mutate(cm_asthma = ifelse(grepl('Asthma', Comorbidities, fixed=TRUE), 1, 0)) %>%
    mutate(cm_copd = ifelse(grepl('Chronic Obstructive Pulmonary Disease (COPD)', Comorbidities, fixed=TRUE), 1, 0)) %>%
    mutate(cm_stroke_any = ifelse((cm_stroke_nos==0) & (cm_ischaemic_stroke==0), 0, 1))

  #save(lab_full, file='data/val.rda')
  lab_full
}

coh.parse_contacts <- function(lab_full,
                               interventions = c('Interventions - Dietetics', 'Interventions - Occupational Therapy',
                                                 'Interventions - Physiotherapy', 'Interventions - Speech and Language Therapy'),
                               int_codes = c('NURSE', 'PT', 'OT', 'SPL', 'DT', 'INFSV', 'CP', 'PC', 'REACT', 'SOCW', 'SUPPW',
                                             'IPMT', 'PT_OT', 'NURSE_PT', 'NURSE_OT', 'OTHER')) {

  print('Linking fictional TRAK Care contact sources...')
  lab_ids <- unique(lab_full$pid)
  trak_pw_test <- tc_pathways %>%
    filter(pid %in% lab_ids) %>%
    mutate(cdt = ymd(dateCreated)) %>%
    mutate(type = 'Inpatient') %>%
    mutate(intervention=careProvider) %>%
    filter(intervention %in% int_codes) %>%
    select(c(pid, cdt, type)) %>%
    distinct()

  trak_ahp_test <- contacts %>%
    filter(pid %in% lab_ids) %>%
    mutate(cdt = ymd(ContactDate)) %>%
    mutate(EpisodeSubType = ifelse((EpisodeSubType == 'Daycase'), 'Inpatient', EpisodeSubType)) %>%
    mutate(EpisodeSubType = ifelse((EpisodeSubType == 'Daycase'), 'Community', EpisodeSubType)) %>%
    filter(EpisodeSubType == 'Inpatient') %>%
    filter(InterventionCategory1 %in% interventions) %>%
    mutate(InterventionCategory1 = str_trim(sub(".*Interventions - ", "", InterventionCategory1))) %>%
    rename(type = EpisodeSubType) %>%
    rename(intervention = InterventionCategory1) %>%
    mutate(intervention=case_when(intervention=='Dietetics'~'DT', intervention=='Physiotherapy'~'PT',
                                  intervention=='Occupational Therapy'~'OT', intervention=='Speech and Language Therapy'~'SPL', TRUE~'OTHER')) %>%
    filter(intervention %in% int_codes) %>%
    select(c(pid, cdt, type)) %>%
    distinct()

  pw_inter <- intersect(trak_pw_test, trak_ahp_test)
  pw_inter <- pw_inter %>%
    drop_na() %>%
    mutate(pr_key = paste(pid, '-', cdt, '-', type))

  trak_pw_test <- trak_pw_test %>%
    drop_na() %>%
    mutate(pr_key = paste(pid, '-', cdt, '-', type))

  trak_ahp_test <- trak_ahp_test %>%
    drop_na() %>%
    mutate(pr_key = paste(pid, '-', cdt, '-', type))

  trak_pw_all <- union(trak_pw_test, trak_ahp_test)
  ahp_ids <- unique(trak_ahp_test$pr_key)
  trakpw_ids <- unique(trak_pw_test$pr_key)
  inter_ids <- unique(pw_inter$pr_key)

  print(sprintf('Parsed: %.0f contacts', dim(trak_pw_all[0]))[1])
  print('Reformatting care provider data..')
  ### Get basic rehab data
  trak_pw_all <- trak_pw_all %>%
    mutate(source = case_when((pr_key %in% inter_ids) ~ 'BOTH', (pr_key %in% trakpw_ids) ~ 'PW_TRAK',
                              (pr_key %in% ahp_ids) ~ 'PW_AHP', TRUE~'Undefined'))

  trak_ahp <- contacts %>%
    filter(pid %in% lab_ids) %>%
    mutate(cdt = ymd(ContactDate)) %>%
    mutate(cduration_ahp = ContactDuration) %>%
    mutate(EpisodeSubType = ifelse((EpisodeSubType == 'Daycase'), 'Inpatient', EpisodeSubType)) %>%
    mutate(EpisodeSubType = ifelse((EpisodeSubType == 'Daycase'), 'Community', EpisodeSubType)) %>%
    filter(EpisodeSubType == 'Inpatient') %>%
    filter(InterventionCategory1 %in% interventions) %>%
    mutate(InterventionCategory1 = str_trim(sub(".*Interventions - ", "", InterventionCategory1))) %>%
    rename(type = EpisodeSubType) %>%
    rename(cintervention_ahp = InterventionCategory1) %>%
    mutate(cintervention_ahp=case_when(cintervention_ahp=='Dietetics'~'DT', cintervention_ahp=='Physiotherapy'~'PT',cintervention_ahp=='Occupational Therapy'~'OT', cintervention_ahp=='Speech and Language Therapy'~'SPL', TRUE~'OTHER')) %>%
    filter(cintervention_ahp %in% int_codes) %>%
    rename(ctime_ahp = ContactTime) %>%
    select(c(pid, cdt, type, ctime_ahp, cduration_ahp, cintervention_ahp, Intervention1, MainReasonActivity, PatientLocationCategory)) %>%
    rename(intervention_details=Intervention1) %>%
    rename(activity=MainReasonActivity) %>%
    rename(location_type=PatientLocationCategory) %>%
    drop_na() %>%
    mutate(pr_key = paste(pid, '-', cdt, '-', type)) %>%
    distinct(pr_key, .keep_all = TRUE)

  trak_pw <- tc_pathways %>%
    filter(pid %in% lab_ids) %>%
    mutate(type = 'Inpatient') %>%
    mutate(cintervention_pw=careProvider) %>%
    filter(cintervention_pw %in% int_codes) %>%
    mutate(cdt = ymd(dateCreated)) %>%
    mutate(ctime_pw = timeCreated) %>%
    mutate(activity_pw = 'Undefined activity') %>%
    mutate(intervention_details_pw = 'No information') %>%
    mutate(location_type_pw = 'Wards / OPDs') %>%
    select(c(pid, cdt, type, ctime_pw, cintervention_pw, intervention_details_pw, activity_pw,
             location_type_pw)) %>%
    drop_na() %>%
    mutate(cduration_pw = NA) %>%
    mutate(pr_key = paste(pid, '-', cdt, '-', type)) %>%
    distinct(pr_key, .keep_all = TRUE)

  trak_pw_all <- left_join(trak_pw_all, trak_ahp, by=c('pid', 'cdt', 'type', 'pr_key'))
  trak_pw_all <- left_join(trak_pw_all, trak_pw, by=c('pid', 'cdt', 'type', 'pr_key'))

  trak_pw_all <- trak_pw_all %>%
    mutate(ctime_ahp = ifelse(is.na(ctime_ahp), '23:59', ctime_ahp)) %>%
    mutate(ctime_pw = ifelse(is.na(ctime_pw), '23:59', ctime_pw)) %>%
    rowwise() %>%
    mutate(ctime = min(ctime_ahp,ctime_pw)) %>%
    mutate(cintervention = case_when((ctime_pw < ctime_ahp)&(!is.na(cintervention_pw)) ~ cintervention_pw,
                                     (ctime_ahp < ctime_pw)&(!is.na(cintervention_ahp)) ~ cintervention_ahp,
                                     (is.na(cintervention_ahp))&(!is.na(cintervention_pw)) ~ cintervention_pw,
                                     (is.na(cintervention_pw))&(!is.na(cintervention_ahp)) ~ cintervention_ahp)) %>%
    ungroup() %>%
    mutate(cintervention = ifelse((is.na(cintervention))&(!is.na(cintervention_pw)), cintervention_pw, cintervention)) %>%
    mutate(cintervention = ifelse((is.na(cintervention))&(!is.na(cintervention_ahp)), cintervention_ahp, cintervention)) %>%
    filter(!is.na(cintervention)) %>%
    mutate(cduration = cduration_ahp) %>%
    mutate(intervention_details = ifelse(is.na(intervention_details), 'No information', intervention_details)) %>%
    mutate(activity = ifelse(is.na(activity), cintervention_pw, 'Undefined activity')) %>%
    mutate(location_type = ifelse(is.na(location_type), 'Wards / OPDs', location_type)) %>%
    select(c(pid, cdt, type, source, ctime, cduration, cintervention, activity, intervention_details, activity, location_type)) %>%
    mutate(cdtt = ymd_hm(paste(cdt, ctime)))

  print(sprintf('Parsed: %.0f contacts', dim(trak_pw_all[0]))[1])
  print('Setting 90-day follow-up and out-of-hours visits...')
  lab_contacts <- left_join(lab_full[, c('pid', 'AdmissionDate', 'DischargeDate',
                                           'AdmissionDateTime', 'DischargeDateTime')], trak_pw_all, by=c('pid'))
  lab_contacts$AdmissionDate <- lubridate::ymd(lab_contacts$AdmissionDate)
  lab_contacts$DischargeDate <- lubridate::ymd(lab_contacts$DischargeDate)
  lab_contacts$DischargeDateTime <- lubridate::ymd_hms(lab_contacts$DischargeDateTime)
  lab_contacts$AdmissionDateTime <- lubridate::ymd_hms(lab_contacts$AdmissionDateTime)
  lab_contacts$cdt <- lubridate::ymd(lab_contacts$cdt)
  lab_contacts <- lab_contacts %>%
    mutate(adm_90d = AdmissionDate + days(90)) %>%
    filter(AdmissionDate < cdt) %>%
    ## Define out of hours rehab
    mutate(out_of_hours = ifelse((hm(as.character(ctime)) > hm('19:00')) | (hm(as.character(ctime)) < hm('07:00')), 1, 0)) %>%
    filter(type=='Inpatient') %>%
    ### Add 90-day contact follow-up
    filter(ymd(cdt) < adm_90d) %>%
    mutate(isFirstWave = ifelse((cdt >= as.Date("2020-03-01")) & (cdt <= as.Date("2020-05-31")), 1, 0)) %>%
    mutate(isSecondWave = ifelse((cdt >= as.Date("2020-09-01")) & (cdt <= as.Date("2021-03-31")), 1, 0)) %>%
    mutate(cov_wave=ifelse(isFirstWave==1, 'First', 'Second'))

  print(sprintf('Parsed: %.0f contacts', dim(lab_contacts[0]))[1])


  lab_contacts <- left_join(lab_contacts, lab_full[, c('pid', 'age_at_first_test', 'sex', 'simd_quint')],
                            by=c('pid'))

  print('Getting flag detailing consecutive contacts within 7 days..')
  ### Filter 2+ 7d contacts
  interv_sel <- lab_contacts %>%
    arrange(pid, cdt) %>%
    mutate(time_to_therapy = as.numeric(difftime(cdt, AdmissionDate, units=c('days'))))

  first_interv <- interv_sel %>%
    arrange(pid, AdmissionDate) %>%
    select(c(pid, time_to_therapy)) %>%
    group_by(pid) %>%
    slice(1)

  second_interv <- interv_sel %>%
    arrange(pid, AdmissionDate) %>%
    select(c(pid, time_to_therapy)) %>%
    group_by(pid) %>%
    rename(time_to_therapy_2=time_to_therapy) %>%
    slice(2)

  lab_contacts <- left_join(lab_contacts, first_interv, by='pid')
  lab_contacts <- left_join(lab_contacts, second_interv, by='pid')
  lab_contacts$time_to_therapy[is.na(lab_contacts$time_to_therapy)] <- -1
  lab_contacts$time_to_therapy_2[is.na(lab_contacts$time_to_therapy_2)] <- -1
  lab_contacts$time_to_therapy <- round(lab_contacts$time_to_therapy)
  lab_contacts$time_to_therapy_2 <- round(lab_contacts$time_to_therapy_2)

  lab_contacts <- lab_contacts %>%
    rowwise() %>%
    mutate(ttt_diff = time_to_therapy_2 - time_to_therapy) %>%
    ungroup() %>%
    mutate(second_contact_7d = ifelse((ttt_diff <= 7)&(time_to_therapy_2!=-1)&(time_to_therapy!=-1), 1, 0)) %>%
    filter(second_contact_7d == 1)

  ## Impute missing contacts
  imp_set <- zoo(lab_contacts$cduration, lab_contacts$cdt)
  lab_contacts <- lab_contacts %>%
    mutate(cduration = na.interpolation(imp_set, option = "linear")) %>%
    mutate(cduration = as.integer(cduration))

  contact_ids <- c(unique(lab_contacts$pid))

  print(sprintf('Parsed: %.0f contacts', dim(lab_contacts[0]))[1])
  print('Generating rehab summary features..')
  ### Create trak pathway rehab features
  interv_sel <- lab_contacts %>%
    arrange(pid, cdt) %>%
    mutate(time_to_therapy = as.numeric(difftime(cdt, AdmissionDate, units=c('days'))))

  first_interv <- interv_sel %>%
    arrange(pid, AdmissionDate) %>%
    select(c(pid, time_to_therapy)) %>%
    group_by(pid) %>%
    slice(1)

  second_interv <- interv_sel %>%
    arrange(pid, AdmissionDate) %>%
    select(c(pid, time_to_therapy)) %>%
    group_by(pid) %>%
    rename(time_to_therapy_2=time_to_therapy) %>%
    slice(2)

  first_type <- interv_sel %>%
    arrange(pid, AdmissionDate) %>%
    select(c(pid, type)) %>%
    group_by(pid) %>%
    slice(1)

  second_type <- interv_sel %>%
    arrange(pid, AdmissionDate) %>%
    select(c(pid, type)) %>%
    group_by(pid) %>%
    rename(type_2=type) %>%
    slice(2)

  first_interv_name <- interv_sel %>%
    arrange(pid, AdmissionDate) %>%
    select(c(pid, cintervention)) %>%
    group_by(pid) %>%
    slice(1)

  second_interv_name <- interv_sel %>%
    arrange(pid, AdmissionDate) %>%
    select(c(pid, cintervention)) %>%
    group_by(pid) %>%
    rename(cintervention_2=cintervention) %>%
    slice(2)

  n_contact_count <- interv_sel %>%
    select(pid) %>%
    group_by(pid) %>%
    summarise(total_count_rehab = n())

  n_ooh_contact_count <- interv_sel %>%
    filter(out_of_hours==1) %>%
    select(pid) %>%
    group_by(pid) %>%
    summarise(total_count_ooh_rehab = n())

  n_contact_minutes <- interv_sel %>%
    select(c(pid, cduration)) %>%
    group_by(pid) %>%
    summarise(total_mins_rehab = sum(cduration))

  n_ooh_contact_minutes <- interv_sel %>%
    filter(out_of_hours==1) %>%
    select(c(pid, cduration)) %>%
    group_by(pid) %>%
    summarise(total_sum_ooh_rehab = sum(cduration))

  ### List of intervention types
  interv_types <- within(interv_sel, {cintervention <- as.character(cintervention);
  rehab_interventions <- ave(cintervention, pid, FUN=toString)})

  interv_types <- interv_types %>%
    select(c(pid, rehab_interventions)) %>%
    group_by(pid) %>%
    slice(1)

  ### Merge all data
  lab_full <- left_join(lab_full, interv_types, by='pid')
  lab_full <- left_join(lab_full, first_type, by='pid')
  lab_full <- left_join(lab_full, second_type, by='pid')
  lab_full <- left_join(lab_full, n_contact_count, by='pid')
  lab_full <- left_join(lab_full, n_ooh_contact_count, by='pid')
  lab_full <- left_join(lab_full, n_contact_minutes, by='pid')
  lab_full <- left_join(lab_full, n_ooh_contact_minutes, by='pid')
  lab_full <- left_join(lab_full, first_interv_name, by='pid')
  lab_full <- left_join(lab_full, second_interv_name, by='pid')
  lab_full <- left_join(lab_full, first_interv, by='pid')
  lab_full <- left_join(lab_full, second_interv, by='pid')
  lab_full$time_to_therapy[is.na(lab_full$time_to_therapy)] <- -1
  lab_full$time_to_therapy_2[is.na(lab_full$time_to_therapy_2)] <- -1
  lab_full$time_to_therapy <- round(lab_full$time_to_therapy)
  lab_full$time_to_therapy_2 <- round(lab_full$time_to_therapy_2)

  lab_full$total_count_rehab[is.na(lab_full$total_count_rehab)] <- 0
  lab_full$total_count_ooh_rehab[is.na(lab_full$total_count_ooh_rehab)] <- 0
  lab_full$cintervention[is.na(lab_full$cintervention)] <- 'No intervention'
  lab_full$cintervention_2[is.na(lab_full$cintervention_2)] <- 'No intervention'
  lab_full$total_sum_ooh_rehab[is.na(lab_full$total_sum_ooh_rehab)] <- 0
  lab_full$total_mins_rehab[is.na(lab_full$total_mins_rehab)] <- 0
  lab_full <- lab_full %>%
    filter(pid %in% contact_ids)

  lab_full <- lab_full %>%
    mutate(n_physio = str_count(lab_full$rehab_interventions, "PT")) %>%
    mutate(n_physio = ifelse(is.na(n_physio), 0, n_physio)) %>%
    mutate(n_occ = str_count(lab_full$rehab_interventions, "OT")) %>%
    mutate(n_occ = ifelse(is.na(n_occ), 0, n_occ)) %>%
    mutate(n_diet = str_count(lab_full$rehab_interventions, "DT")) %>%
    mutate(n_diet = ifelse(is.na(n_diet), 0, n_diet)) %>%
    mutate(n_spl = str_count(lab_full$rehab_interventions, "SPL")) %>%
    mutate(n_spl = ifelse(is.na(n_spl), 0, n_spl)) %>%
    mutate(n_nurse = str_count(lab_full$rehab_interventions, "NURSE")) %>%
    mutate(n_nurse = ifelse(is.na(n_nurse), 0, n_nurse))

  ### Add admission start, end flags and deaths to contacts data
  print('Appending additional events to contacts data..')
  lab_contacts <- left_join(lab_contacts, deaths, by='pid')

  lab_adm <- lab_contacts %>%
    arrange(pid, AdmissionDateTime) %>%
    group_by(pid) %>%
    slice(1) %>%
    mutate(intervention_details='Point of admission') %>%
    mutate(cintervention='ADM-START') %>%
    mutate(cduration=0) %>%
    mutate(cdtt=AdmissionDateTime) %>%
    mutate(cdt=AdmissionDate)

  lab_disch <- lab_contacts %>%
    arrange(pid, DischargeDateTime) %>%
    group_by(pid) %>%
    slice(1) %>%
    mutate(intervention_details='Point of discharge') %>%
    mutate(cintervention='DISCH-END') %>%
    mutate(cduration=0) %>%
    mutate(cdtt=DischargeDateTime) %>%
    mutate(cdt=DischargeDate) %>%
    filter(ymd(cdt) < adm_90d)

  lab_death <- lab_contacts %>%
    filter(!is.na(date_of_death)) %>%
    filter(ymd(date_of_death) <= ymd(adm_90d)) %>%
    arrange(pid, date_of_death) %>%
    group_by(pid) %>%
    slice(1) %>%
    mutate(intervention_details='Point of death') %>%
    mutate(cintervention='DEATH') %>%
    mutate(cduration=0) %>%
    mutate(cdtt=ymd_hms(paste(date_of_death, "00:00:00", sep=" "))) %>%
    mutate(cdt=date_of_death)

  lab_contacts_full <- rbind(lab_contacts, lab_adm)
  lab_contacts_full <- rbind(lab_contacts_full, lab_disch)
  lab_contacts_full <- rbind(lab_contacts_full, lab_death)

  print('Reformatting contacts data and cohorts file..')
  ### Omit missing dates
  lab_contacts_full <- lab_contacts_full %>%
    mutate(cdtt_new = parse_date_time(cdtt, orders = c("ymd HMS"))) %>%
    mutate(cdtt_new = as.character(cdtt_new)) %>%
    mutate(cdtt_new = ifelse(is.na(cdtt_new), paste0(cdtt, "00:00:00"), cdtt_new)) %>%
    mutate(cdtt=cdtt_new) %>%
    select(-cdtt_new) %>%
    filter(!is.na(AdmissionDateTime)) %>%
    filter(!is.na(DischargeDateTime)) %>%
    filter(!is.na(cdt)) %>%
    filter(!is.na(cdtt)) %>%
    filter(!is.na(adm_90d))

  ### Set comorbidity flags
  lab_full <- lab_full %>%
    mutate(cm_chd=ifelse(cm_chd==0, 'N', 'Y')) %>%
    mutate(cm_stroke_any=ifelse(cm_stroke_any==0, 'N', 'Y')) %>%
    mutate(cm_mi=ifelse(cm_mi==0, 'N', 'Y')) %>%
    mutate(cm_diabetes=ifelse(cm_diabetes==0, 'N', 'Y')) %>%
    mutate(cm_obesity=ifelse(cm_obesity==0, 'N', 'Y')) %>%
    mutate(cm_hip_fracture=ifelse(cm_hip_fracture==0, 'N', 'Y')) %>%
    mutate(cm_dementia=ifelse(cm_dementia==0, 'N', 'Y')) %>%
    mutate(cm_delirium=ifelse(cm_delirium==0, 'N', 'Y')) %>%
    mutate(cm_depression=ifelse(cm_depression==0, 'N', 'Y')) %>%
    mutate(cm_asthma=ifelse(cm_asthma==0, 'N', 'Y')) %>%
    mutate(cm_copd=ifelse(cm_copd==0, 'N', 'Y')) %>%
    mutate(isFirstWave = ifelse((AdmissionDate >= as.Date("2020-03-01")) & (AdmissionDate <= as.Date("2020-05-31")), 1, 0)) %>%
    mutate(isSecondWave = ifelse((AdmissionDate >= as.Date("2020-09-01")) & (AdmissionDate <= as.Date("2021-03-31")), 1, 0)) %>%
    mutate(cov_wave=ifelse(isFirstWave==1, 'First', 'Second'))


  ### Set Length of stay
  lab_full <- lab_full %>%
    mutate(index_LOS = as.integer(difftime(DischargeDate, AdmissionDate, units='days')))


  ### Set Rehab to Length-of-stay ratio
  lab_full <- lab_full %>%
    rowwise() %>%
    mutate(rehab_to_los = round(total_mins_rehab / round(index_LOS, 0), 2)) %>%
    mutate(rehab_to_los = ifelse(is.na(rehab_to_los), 0, rehab_to_los)) %>%
    mutate(rehab_to_los = ifelse(is.infinite(rehab_to_los), total_mins_rehab, rehab_to_los)) %>%
    mutate(total_sum_ooh_rehab = as.numeric(total_sum_ooh_rehab))



  ### Add age groups and quintiles for rehab minutes
  lab_contacts_full <- lab_contacts_full %>%
    mutate(age_group=case_when(age_at_first_test<65~'18-64', (age_at_first_test>=65 & age_at_first_test<80)~'65-79',
                               age_at_first_test>79~'80+'))


  lab_full$rmin_group <- as.numeric(cut2(lab_full$total_mins_rehab, g=5))
  lab_contacts_full <- left_join(lab_contacts_full, lab_full[,c('pid', 'rmin_group')], by='pid')

  print('Complete..')
  save(lab_full, file='data/cov_coh.rda')
  save(lab_contacts_full, file='data/cov_ct.rda')

  print(sprintf('Complete, parsed: %.0f contacts for %.0f patients.', dim(lab_contacts_full[0]), dim(lab_full[0]))[1])
}

coh.adjust_distribution <- function(adjust_days=300, adj_prop=1.1, adj_los=200, adj_dlos=250) {
  load(file='data/cov_coh.rda')
  load(file='data/cov_ct.rda')

  print('Adjusting randomised episodes to balance population..')
  lab_adj <- lab_full %>%
    slice_sample(n=round(dim(lab_full)[1]/adj_prop)) %>%
    mutate(pid = pid+100000) %>%
    mutate(AdmissionDate = as.Date(AdmissionDate) + days(adjust_days)) %>%
    #mutate(AdmissionDateTime = list(rtime(round(dim(lab_full)[1]/adj_prop)))) %>%
    #mutate(AdmissionDateTime = ymd_hm(paste(AdmissionDate, AdmissionDateTime))) %>%
    mutate(collectionDate = as.Date(collectionDate) + days(adjust_days)) %>%
    filter(AdmissionDate <= DischargeDate)

  lab_mg <- rbind(lab_adj, lab_full)
  lab_mg <- lab_mg %>%
    select(-c(DischargeDateTime, AdmissionDateTime)) %>%
    select(-c(isFirstWave, isSecondWave)) %>%
    mutate(isFirstWave = ifelse((AdmissionDate >= as.Date("2020-03-01")) & (AdmissionDate <= as.Date("2020-05-31")), 1, 0)) %>%
    mutate(isSecondWave = ifelse((AdmissionDate >= as.Date("2020-09-01")) & (AdmissionDate <= as.Date("2021-03-31")), 1, 0)) %>%
    mutate(cov_wave = ifelse(isFirstWave == 1, 'First', 'None')) %>%
    mutate(cov_wave = ifelse(isSecondWave == 1, 'Second', cov_wave)) %>%
    mutate(DischargeDate = ifelse((index_LOS>adj_los)&(isFirstWave==1), (as.Date(DischargeDate) - days(adj_dlos)), DischargeDate)) %>%
    mutate(DischargeDate = as.Date(DischargeDate)) %>%
    filter(AdmissionDate <= DischargeDate) %>%
    mutate(index_LOS = as.integer(difftime(DischargeDate, AdmissionDate, units='days'))) %>%
    mutate(total_mins_rehab = round(total_mins_rehab / 10)) %>%
    mutate(rehab_to_los = round(total_mins_rehab / round(index_LOS, 0), 2)) %>%
    mutate(rehab_to_los = ifelse(is.na(rehab_to_los), 0, rehab_to_los)) %>%
    mutate(rehab_to_los = ifelse(is.infinite(rehab_to_los), total_mins_rehab, rehab_to_los))

  print('Adjusting contacts to balance population..')
  lab_cts_adj <- lab_contacts_full %>%
    slice_sample(n=round(dim(lab_contacts_full)[1]/adj_prop)) %>%
    mutate(pid = pid+100000) %>%
    mutate(cdt = as.Date(cdt) + days(adjust_days))

  lab_contacts_fn <- rbind(lab_cts_adj, lab_contacts_full)

  lab_contacts_fn <- lab_contacts_fn %>%
    mutate(ctime = rtime(round(dim(lab_contacts_fn)[1]))) %>%
    select(-cdtt) %>%
    mutate(cdtt = ymd_hm(paste(cdt, ctime))) %>%
    mutate(wv1_contact = ifelse((cdt >= as.Date("2020-03-01")) & (cdt <= as.Date("2020-05-31")), 1, 0)) %>%
    mutate(wv2_contact = ifelse((cdt >= as.Date("2020-09-01")) & (cdt <= as.Date("2021-03-31")), 1, 0)) %>%
    mutate(cov_wave = ifelse(wv1_contact == 1, 'First', 'None')) %>%
    mutate(cov_wave = ifelse(wv2_contact == 1, 'Second', cov_wave))

  print('Adjustment complete..')
  save(lab_mg, file='data/cov_coh.rda')
  save(lab_contacts_fn, file='data/cov_ct.rda')
}
