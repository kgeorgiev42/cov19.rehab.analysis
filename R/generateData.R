#' Generates dummy data for COVID-19 rehabilitation analysis.
#'
#'

gen_cov_data <- function() {
  set.seed(42)
  print('Generating dummy data for analysis..')
  print('Preparing COVID-19 episodes...')
  episodes <- gen_cov_data.cov19_episodes()
  print('Preparing COVID-19 lab tests...')
  lab_tests <- gen_cov_data.cov19_lab_tests()
  print('Preparing AHP contacts...')
  contacts <- gen_cov_data.contacts()
  print('Preparing Care pathway data...')
  trak_care <- gen_cov_data.trak_care_pathways()
  print('Preparing Phenotype data...')
  smr_phen <- gen_cov_data.smr_phenotypes()
  print('Preparing Admission episodes...')
  admissions <- gen_cov_data.trak_admissions()
  print('Preparing Discharge episodes...')
  discharges <- gen_cov_data.trak_discharges()
  print('Preparing Deaths data...')
  deaths <- gen_cov_data.deaths()
}

gen_cov_data.cov19_episodes <- function(nrows=200000, fm_prob=0.54, simd_probs=c(0.14,0.24,0.17,0.18,0.26),
                           age_mp1=c(70,7), age_mp2=c(85,8)) {
   age_m1 <- rnorm(nrows/2, age_mp1[1], age_mp1[2])
   age_m2 <- rnorm(nrows/2, age_mp2[1], age_mp2[2])
   ages <- c(age_m1, age_m2)
   episodes <- data.frame(pid=seq(1,nrows),
                          sex=rbinom(nrows,1,fm_prob),
                          simd_quint=sample(c(1,2,3,4,5), size=nrows, replace=TRUE, prob = simd_probs),
                          age_at_first_test=ages)

   episodes <- episodes %>%
     arrange(pid) %>%
     mutate(sex=ifelse(sex==0, 'F', 'M')) %>%
     mutate(age_at_first_test=as.integer(age_at_first_test))

   save(episodes, file = "data/episodes.rda")

}

gen_cov_data.cov19_lab_tests <- function(nrows=200000, n_pts=100000, max_tppt=4, ds_prob=0.51,
                                    tp_prob=1.0, start='2020-03-01', end='2021-03-31',
                                    follow_up=365) {

  # Create test IDs linked to each patient
  test_ids <- lapply(1:n_pts, function(pt) {
    num_tests <- sample(1:max_tppt, 1)
    tests <- sample(1:nrows, num_tests, replace = FALSE)
    data.frame(pid = rep(pt, num_tests), lab_id = tests)
  })
  lab_tests <- do.call(rbind, test_ids)
  lab_tests <- lab_tests %>%
    mutate(dataSource=rbinom(dim(lab_tests)[1], 1, ds_prob)) %>%
    mutate(dataSource=ifelse(dataSource==0, 'Test by Lighthouse Laboratory', 'iLabs')) %>%
    mutate(testPositive=rbinom(dim(lab_tests)[1], 1, tp_prob)) %>%
    mutate(testPositive=ifelse(testPositive==0, 'N', 'Y')) %>%
    mutate(collectionDate=rdate(dim(lab_tests)[1], min=start, max=end)) %>%
    mutate(collectionDate=as.Date(sapply(collectionDate, rdate_increase, end=follow_up))) %>%
    filter(collectionDate < as.Date(end)) %>%
    arrange(pid, collectionDate)

  save(lab_tests, file = "data/lab_tests.rda")

}

gen_cov_data.contacts <- function(n_pts=25000, max_cppt=60,
                                  ip_prob=0.6, int_probs=c(0.63,0.20,0.05,0.12),
                                  loc_probs=c(0.55, 0.20, 0.09, 0.16),
                                  duration_params=c(65,350,5,0.09),
                                  dgr_probs=c(0.96, 0.04),
                                  act_probs=c(0.01, 0.002, 0.31, 0.01, 0.07, 0.01,
                                              0.02, 0.003, 0.11, 0.03, 0.12, 0.001,
                                              0.02, 0.007, 0.007, 0.004, 0.18, 0.004,
                                              0.004, 0.04, 0.02),
                                  start='2020-03-01', end='2021-03-31', follow_up=365) {

  # Create AHP contacts linked to each patient
  contact_ids <- lapply(1:n_pts, function(pt) {
    num_cts <- sample(1:max_cppt, 1)
    data.frame(pid = rep(pt, num_cts))
  })
  contacts <- do.call(rbind, contact_ids)
  ct_durations <- rpois(dim(contacts)[1], lambda = duration_params[1]) - duration_params[3]
  ct_durations <- pmin(pmax(ct_durations, 0), duration_params[2])

  contacts <- contacts %>%
    mutate(EpisodeSubType=rbinom(dim(contacts)[1],1, ip_prob)) %>%
    mutate(EpisodeSubType=ifelse(EpisodeSubType==0, 'Community', 'Inpatient')) %>%
    mutate(MainReasonActivity=sample(c('Physiotherapy intervention','Occupational Therapy intervention',
                                       'Speech and Language intervention', 'Dietetic intervention'),
                                     size=dim(contacts)[1], replace=TRUE, prob=int_probs)) %>%
    mutate(AHPProfession=case_when(MainReasonActivity=='Physiotherapy intervention'~'Physiotherapists',
                                   MainReasonActivity=='Occupational Therapy intervention'~'Occupational therapists',
                                   MainReasonActivity=='Speech and Language intervention'~'Speech and Language therapists',
                                   MainReasonActivity=='Dietetic intervention'~'Dieticians')) %>%
    mutate(PatientLocationCategory=sample(c('Wards / OPDs','NHS hospitals (inc day)',
                                            'Health Centres / GP Surgeries', 'Other Locations (generic)'),
                                          size=dim(contacts)[1], replace=TRUE, prob=loc_probs)) %>%
    mutate(ContactDuration=ct_durations) %>%
    mutate(ContactDate=rdate(dim(contacts)[1], min=start, max=end)) %>%
    mutate(ContactTime=rtime(dim(contacts)[1])) %>%
    mutate(ContactTime=ifelse(ContactTime=='NA:NA', '00:00', ContactTime)) %>%
    mutate(ContactDayGroup=sample(c('Mon-Fri', 'Sat-Sun'), size=dim(contacts)[1], replace=TRUE, prob=dgr_probs)) %>%
    mutate(InterventionCategory1=case_when(MainReasonActivity=='Physiotherapy intervention'~'Interventions - Physiotherapy',
                                           MainReasonActivity=='Occupational Therapy intervention'~'Interventions - Occupational Therapy',
                                           MainReasonActivity=='Speech and Language intervention'~'Interventions - Speech and Language Therapy',
                                           MainReasonActivity=='Dietetic intervention'~'Interventions - Dietetics')) %>%
    mutate(Intervention1=sample(c('A-DIET-NUTR', 'A-PHYS-FUNC', 'A-PHYS-GEN', 'A-PHYS-MOB',
                                  'A-SPL-EDS', 'A-SPL-SLC', 'ADMIN-REFTRT', 'ADV-GEN',
                                  'T-DIET-ETF', 'T-DIET-FF', 'T-DIET-ONS', 'T-DIET-TM',
                                  'T-DIET-TPN', 'T-PHYS-EXER', 'T-PHYS-GB', 'T-PHYS-MOB',
                                  'T-PHYS-RESP', 'T-PHYS-STAIRS', 'T-PHYS-TRANSFER', 'T-SPL-EDS',
                                  'T-SPL-SLC'), size=dim(contacts)[1], replace=TRUE, prob=act_probs)) %>%
    mutate(ContactDate=as.Date(sapply(ContactDate, rdate_increase, end=follow_up))) %>%
    arrange(pid, ContactDate, ContactTime) %>%
    filter(ContactDate < as.Date(end))

    ### Add NAs
    na_indices <- sample(1:dim(contacts)[1], round(dim(contacts)[1] * duration_params[4]))
    contacts[na_indices, 'ContactDuration'] <- NA

  save(contacts, file = "data/contacts.rda")

}

gen_cov_data.trak_care_pathways <- function(n_pts=25000, max_cppt=60,
                                              cp_probs=c(0.003, 0.08, 0.004, 0.0001,
                                                         0.26, 0.03, 0.0001, 0.59, 0.03),
                                            start='2020-03-01', end='2021-03-31',
                                            follow_up=365) {
  # Create TrakCare contacts linked to each patient
  contact_ids <- lapply(1:n_pts, function(pt) {
    num_cts <- sample(1:max_cppt, 1)
    data.frame(pid = rep(pt, num_cts))
  })
  tc_pathways <- do.call(rbind, contact_ids)

  tc_pathways <- tc_pathways %>%
    mutate(dateCreated=rdate(dim(tc_pathways)[1], min=start, max=end)) %>%
    mutate(timeCreated=rtime(dim(tc_pathways)[1])) %>%
    mutate(careProvider=sample(c('CP', 'DT', 'INFSV', 'IPMT',
                                  'NURSE', 'OT', 'PC', 'PT',
                                  'SPL'), size=dim(tc_pathways)[1], replace=TRUE, prob=cp_probs)) %>%
    mutate(dateCreated=as.Date(sapply(dateCreated, rdate_increase, end=follow_up))) %>%
    filter(dateCreated < as.Date(end)) %>%
    arrange(pid, dateCreated, timeCreated)

  save(tc_pathways, file = "data/trak_care.rda")

}

gen_cov_data.smr_phenotypes <- function(n_pts=42000, max_cppt=12,
                                            cp_probs=c(0.013, 0.011, 0.016, 0.01,
                                                       0.006, 0.008, 0.015, 0.007,
                                                       0.006, 0.008, 0.9),
                                        start='1996-01-01', end='2021-03-31') {
  # The Scottish Morbidity Records datasets typically contain ICD-10 or other coded long-term conditions
  # This is a simplified dataset that uses a phenotype structure based on the HDRUK Phenotype Library
  # Create fictional SMR coded events
  smr_ids <- lapply(1:n_pts, function(pt) {
    num_cond <- sample(1:max_cppt, 1)
    data.frame(pid = rep(pt, num_cond))
  })
  smr_phenotypes <- do.call(rbind, smr_ids)
  smr_phenotypes <- smr_phenotypes %>%
    mutate(EventDate=rdate(dim(smr_phenotypes)[1], min=start, max=end)) %>%
    mutate(EventEndDate=as.Date(sapply(EventDate, rdate_increase))) %>%
    mutate(PhenotypeName=sample(c('Asthma', 'COPD', 'Coronary Heart Disease', 'Delirium',
                                 'Dementia', 'Depression', 'Diabetes', 'Stroke', 'Obesity',
                                 'Myocardial Infarction', 'None'), size=dim(smr_phenotypes)[1], replace=TRUE, prob=cp_probs)) %>%
    filter(PhenotypeName!='None')

  save(smr_phenotypes, file = "data/smr_phenotype.rda")
}

gen_cov_data.trak_admissions <- function(nrows=120000, n_pts=100000, max_admppt=4,
                                         start_wv1='2020-03-01', end_wv1='2020-05-31',
                                         start_wv2='2020-09-01', end_wv2='2021-03-31',
                                         follow_up=365) {
  # Create fictional admission IDs
  adm_ids <- lapply(1:n_pts, function(pt) {
    num_adm <- sample(1:max_admppt, 1)
    admissions <- sample(1:nrows, num_adm, replace = FALSE)
    data.frame(pid = rep(pt, num_adm), episode = admissions)
  })
  admissions <- do.call(rbind, adm_ids)
  admissions_wv1 <- admissions %>%
    mutate(AdmissionDate=rdate(dim(admissions)[1], min_date=start_wv1, max_date=end_wv1)) %>%
    mutate(AdmissionTime=rtime(dim(admissions)[1])) %>%
    mutate(AdmissionDate=as.Date(sapply(AdmissionDate, rdate_increase, end=follow_up))) %>%
    arrange(pid, AdmissionDate) %>%
    distinct(pid, AdmissionDate, .keep_all=TRUE) %>%
    filter(AdmissionDate <= as.Date(end_wv1))

  admissions_wv2 <- admissions %>%
    mutate(AdmissionDate=rdate(dim(admissions)[1], min_date=start_wv2, max_date=end_wv2)) %>%
    mutate(AdmissionTime=rtime(dim(admissions)[1])) %>%
    mutate(AdmissionDate=as.Date(sapply(AdmissionDate, rdate_increase, end=follow_up))) %>%
    arrange(pid, AdmissionDate) %>%
    distinct(pid, AdmissionDate, .keep_all=TRUE) %>%
    filter(AdmissionDate <= as.Date(end_wv2))

  admissions <- rbind(admissions_wv1, admissions_wv2)
  save(admissions, file = "data/trak_adm.rda")

}

gen_cov_data.trak_discharges <- function(end_date='2021-03-31', follow_up=30) {

  load(file='data/trak_adm.rda')
  discharges <- data.frame(admissions)

  discharges <- discharges %>%
    mutate(DischargeDate=as.Date(sapply(AdmissionDate, rdate_increase, end=follow_up))) %>%
    mutate(DischargeTime=rtime(dim(discharges)[1])) %>%
    mutate(DischargeDate=ifelse(DischargeDate > end_date, NA, DischargeDate)) %>%
    filter(!is.na(DischargeDate)) %>%
    mutate(DischargeDate=as.Date(DischargeDate)) %>%
    mutate(DischargeDate = ifelse((lead(AdmissionDate) < DischargeDate)&(lead(pid)==pid), as.Date(lead(AdmissionDate)) - days(1),
                                  DischargeDate)) %>%
    mutate(DischargeDate=as.Date(DischargeDate)) %>%
    filter(!is.na(DischargeDate)) %>%
    select(-c(AdmissionDate, AdmissionTime))

  save(discharges, file = "data/trak_disch.rda")

}

gen_cov_data.deaths <- function(nrows=50000, end_date='2021-03-31', follow_up=365) {

  load(file='data/trak_disch.rda')
  load(file='data/smr_phenotype.rda')
  load(file='data/trak_care.rda')
  load(file='data/contacts.rda')
  load(file='data/episodes.rda')
  load(file='data/lab_tests.rda')

  # Create fictional death events
  deaths <- discharges %>%
    arrange(pid, desc(DischargeDate)) %>%
    group_by(pid) %>%
    slice(1) %>%
    ungroup() %>%
    slice_sample(n = nrows) %>%
    mutate(date_of_death=as.Date(sapply(DischargeDate, rdate_increase, end=follow_up))) %>%
    mutate(date_of_death=ifelse(date_of_death>end_date, NA, date_of_death)) %>%
    filter(!is.na(date_of_death)) %>%
    mutate(date_of_death=as.Date(date_of_death)) %>%
    select(-c(DischargeDate, DischargeTime, episode))

  smr_pheno <- smr_phenotypes %>%
    arrange(pid, desc(EventEndDate)) %>%
    group_by(pid) %>%
    slice(1) %>%
    ungroup() %>%
    select(c('pid', 'EventDate', 'EventEndDate'))
  deaths <- left_join(deaths, smr_pheno, by='pid')
  tc_p <- tc_pathways %>%
    arrange(pid, desc(dateCreated)) %>%
    group_by(pid) %>%
    slice(1) %>%
    ungroup() %>%
    select(c('pid', 'dateCreated'))
  deaths <- left_join(deaths, tc_p, by='pid')
  lbs <- lab_tests %>%
    arrange(pid, desc(collectionDate)) %>%
    group_by(pid) %>%
    slice(1) %>%
    ungroup() %>%
    select(c('pid', 'collectionDate'))
  deaths <- left_join(deaths, lbs, by='pid')
  cts <- contacts %>%
    arrange(pid, desc(ContactDate)) %>%
    group_by(pid) %>%
    slice(1) %>%
    ungroup() %>%
    select(c('pid', 'ContactDate'))
  deaths <- left_join(deaths, cts, by='pid')
  deaths <- deaths %>%
    filter(date_of_death >= as.Date(EventEndDate)) %>%
    filter(date_of_death >= as.Date(dateCreated)) %>%
    filter(date_of_death >= as.Date(collectionDate)) %>%
    filter(date_of_death >= as.Date(ContactDate)) %>%
    select(-c(EventDate, EventEndDate, collectionDate, ContactDate, dateCreated))

  save(deaths, file = "data/deaths.rda")
}



