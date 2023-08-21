pm.get_summary_table <- function(out_path='results/pm_baseline.html', positive_only=1) {
  print('Getting baseline characteristics for COVID-19 positive cohort..')
  load(file='data/cov_coh.rda')
  ## Get positive cases only
  if (positive_only) {
    lab_mg <- lab_mg %>%
      filter(testPositive_10d==1)
  }

  covid_pts_tbl <- lab_mg %>%
    select(c(age_at_first_test, sex, simd_quint, index_LOS, cov_wave, mortality_1yr,
             total_mins_rehab, total_count_rehab, total_count_ooh_rehab, n_nurse, n_physio,
             n_occ, n_diet, n_spl, rehab_to_los, time_to_therapy, time_to_therapy_2)) %>%
    tbl_summary(
      by=cov_wave,
      statistic = list(age_at_first_test~"{median} ({p25}, {p75})", index_LOS~"{median} ({p25}, {p75})",
                       total_mins_rehab~"{median} ({p25}, {p75})", total_count_rehab~"{median} ({p25}, {p75})",
                       total_count_ooh_rehab~"{mean} ({sd})", n_nurse~"{mean} ({sd})", n_physio~"{mean} ({sd})",
                       n_occ~"{mean} ({sd})", n_diet~"{mean} ({sd})",
                       n_spl~"{mean} ({sd})", rehab_to_los~"{median} ({p25}, {p75})",
                       rehab_to_los~"{median} ({p25}, {p75})", time_to_therapy~"{median} ({p25}, {p75})",
                       time_to_therapy_2~"{median} ({p25}, {p75})"),
      type = list(total_count_ooh_rehab ~ 'continuous', n_occ~'continuous', n_spl~'continuous'),
      digits = all_continuous()~2,
      missing="no",
      label=list(age_at_first_test~'Age at first test', sex ~ 'Sex', index_LOS ~ 'Length of stay (days)',
                 simd_quint ~ 'SIMD (Quintiles)',
                 mortality_1yr ~ 'Died (1-year all-cause mortality)',
                 time_to_therapy~'Time to first contact (days)',
                 time_to_therapy_2~'Time to second consecutive contact (days)',
                 total_mins_rehab~'Total minutes of rehabilitation',
                 total_count_rehab~'Total number of contacts',
                 total_count_ooh_rehab~'Total number of out-of-hours contacts',
                 n_nurse~'Nursing',
                 n_physio~'Physiotherapy',
                 n_occ~'Occupational Therapy', n_spl~'Speech and Language Therapy', n_diet~'Dietetics',
                 rehab_to_los~'Minutes of rehabiltation per day of hospitalisation')
    ) %>%
    add_p(pvalue_fun= ~style_pvalue(.x, digits=3)) %>%
    bold_labels() %>%
    modify_spanning_header(c('stat_1', 'stat_2')~"**COVID-19 Wave**")

  print(paste0('Table exported to ', out_path))
  covid_pts_tbl %>% as_gt() %>% gtsave(out_path)
}

pm.generate_event_log <- function(out_path='data/pm_event_log.rda', positive_only=1) {
  load(file='data/cov_ct.rda')
  print('Generating event log data..')

  ## Get positive cases only
  if (positive_only) {
    lab_contacts_fn <- left_join(lab_contacts_fn, lab_mg[,c('pid', 'testPositive_10d')], by='pid')
    lab_contacts_fn <- lab_contacts_fn %>%
      filter(testPositive_10d==1)
  }

  event_log_all <- lab_contacts_fn %>%
    rename(start=cdtt) %>%
    mutate(start=as_datetime(start)) %>%
    mutate(complete=as_datetime(start) + minutes(cduration)) %>%
    mutate(resource=location_type)

  ### Setting event log on the specialist level
  alog_all_spec <- event_log_all %>%
    activitylog(case_id='pid', activity_id='cintervention',
                timestamps=c('start', 'complete'),
                resource_id = 'resource')

  print('Printing summary for Wave 1..')
  alog_all_spec %>%
    filter(wv1_contact==1) %>%
    summary

  print('Printing summary for Wave 2..')
  alog_all_spec %>%
    filter(wv2_contact==1) %>%
    summary

  print(paste0('Event log exported to ', out_path))
  save(alog_all_spec, file=out_path)
}

pm.plot_frequent_traces <- function(in_path='data/pm_event_log.rda',
                                    out_path='results/freq_map.pdf', maf=5, minf=50,
                                    wave='Second') {
  load(file=in_path)
  print('Plotting frequent traces as process map..')
  ### Set COVID-19 wave
  if(wave == 'First'){
    prm <- alog_all_spec %>%
      filter(cov_wave=='First')
  }
  else {
    prm <- alog_all_spec %>%
      filter(cov_wave=='Second')
  }

  prm <- prm %>%
    filter_infrequent_flows(min_n=minf) %>%
    filter_activity_frequency(interval=c(maf, NA)) %>%
    ## Remove straight discharges/admission-to-discharge and vice-versa to make map look more realistic
    filter_trim(start_activities = "ADM-START", end_activities=c("DISCH-END", "SPL", "PT", "OT", "DT", "PC", "CP", "INFSV", "IPMT", "NURSE")) %>%
    filter_precedence(antecedents = c("ADM-START"), consequents = c("DISCH-END"), precedence_type = 'directly_follows', filter_method = 'none') %>%
    process_map(render=F, type_nodes=frequency('absolute-case'), type_edges=frequency('relative-case'),
                sec_edges=performance(median,'days'), rankdir = 'LR')

  prm$nodes_df$label <- str_replace_all(prm$nodes_df$label, c('Start'='PWSTART', 'End'='EOFP'))
  export_graph(prm, paste0(out_path), title='Frequent traces.')
  print(paste0('Process map exported to ', out_path))
}

pm.plot_path_for_age_group <- function(in_path='data/pm_event_log.rda',
                                    out_path='results/freq_map_subage.pdf', maf=5, minf=50,
                                    wave='Second', age_group='80+') {
  ### Default plots pathways in Wave 1 for patients aged 80 or above
  load(file=in_path)
  print('Plotting process map subset..')
  ### Set COVID-19 wave
  if(wave == 'First'){
    prm <- alog_all_spec %>%
      filter(cov_wave=='First')
  }
  else {
    prm <- alog_all_spec %>%
      filter(cov_wave=='Second')
  }

  prm <- prm %>%
    filter(age_group==age_group) %>%
    filter_infrequent_flows(min_n=minf) %>%
    filter_activity_frequency(interval=c(maf, NA)) %>%
    ## Remove straight discharges/admission-to-discharge and vice-versa to make map look more realistic
    filter_trim(start_activities = "ADM-START", end_activities=c("DISCH-END", "SPL", "PT", "OT", "DT", "PC", "CP", "INFSV", "IPMT", "NURSE")) %>%
    filter_precedence(antecedents = c("ADM-START"), consequents = c("DISCH-END"), precedence_type = 'directly_follows', filter_method = 'none') %>%
    process_map(render=F, type_nodes=frequency('absolute-case'), type_edges=frequency('relative-case'),
                sec_edges=performance(median,'days'), rankdir = 'LR')

  prm$nodes_df$label <- str_replace_all(prm$nodes_df$label, c('Start'='PWSTART', 'End'='EOFP'))
  export_graph(prm, paste0(out_path), title='Frequent traces (subset by age group).')
  print(paste0('Process map exported to ', out_path))
}

pm.plot_precedence_matrix <- function(in_path='data/pm_event_log.rda',
                          out_path='results/prec_matrix.png',
                          wave='Second') {
  load(file=in_path)
  print('Plotting precedence matrix..')
  ### Set COVID-19 wave
  if(wave == 'First'){
    prm <- alog_all_spec %>%
      filter(cov_wave=='First')
  }
  else {
    prm <- alog_all_spec %>%
      filter(cov_wave=='Second')
  }

  png(out_path, height=600, width=600)
  prm_p <- prm %>%
    ## Remove straight discharges/admission-to-discharge and vice-versa to make map look more realistic
    filter_trim(start_activities = "ADM-START", end_activities=c("DISCH-END", "SPL", "PT", "OT", "DT", "PC", "CP", "INFSV", "IPMT", "NURSE")) %>%
    filter_precedence(antecedents = c("ADM-START"), consequents = c("DISCH-END"), precedence_type = 'directly_follows', filter_method = 'none') %>%
    to_eventlog() %>%
    precedence_matrix(type='relative')
  print(prm_p %>% plot())
  dev.off()
  print(paste0('Precedence matrix exported to ', out_path))
}

pm.plot_throughput_time <- function(in_path='data/pm_event_log.rda',
                                      out_path='results/throughput_time.png',
                                    int_min=NA, int_max=90, units='days') {
  ### Plots the cumulative pathway completion time over patient groups
  load(file=in_path)
  print('Plotting throughput time by age groups..')

  png(out_path, width=3.25,
      height=3.75,
      units="in",
      res=1200,
      pointsize=4)
  prm_p <- alog_all_spec %>%
    mutate(start = as.Date.POSIXct(start)) %>%
    mutate(complete = as.Date.POSIXct(complete)) %>%
    filter(cov_wave != 'None') %>%
    filter_throughput_time(units=units, interval=c(int_min, int_max)) %>%
    group_by(cov_wave, age_group) %>%
    throughput_time('log', units=units)

  print(prm_p %>% plot())
  dev.off()
  print(paste0('Throughput time plot exported to ', out_path))
}

pm.plot_idle_time <- function(in_path='data/pm_event_log.rda',
                                    out_path='results/idle_time.png',
                                    int_min=NA, int_max=200, units='days') {
  ### Plots the cumulative time of 'rehabilitation-related' inactivity across the pathway
  load(file=in_path)
  print('Plotting idle time by age groups..')

  png(out_path, width=3.25,
      height=3.75,
      units="in",
      res=1200,
      pointsize=4)
  prm_p <- alog_all_spec %>%
    mutate(start = as.Date.POSIXct(start)) %>%
    mutate(complete = as.Date.POSIXct(complete)) %>%
    filter(cov_wave != 'None') %>%
    group_by(cov_wave, age_group) %>%
    filter_processing_time(units=units, interval=c(int_min, int_max)) %>%
    idle_time('log', units=units)

  print(prm_p %>% plot())
  dev.off()
  print(paste0('Idle time plot exported to ', out_path))
}

pm.plot_dotted_chart <- function(in_path='data/pm_event_log.rda',
                                 out_path='results/dotted_chart.pdf') {
  ### Plots a cumulative dotted chart of patient events stratified by COVID-19 wave
  load(file=in_path)
  print('Plotting event log dotted chart..')

  alog_plot <- alog_all_spec %>%
    filter(cov_wave != 'None') %>%
    select(-c(activity)) %>%
    mutate(start = as.Date(start)) %>%
    ### Correct adjusted pids from randomisation
    mutate(pid = ifelse((pid>=100000)&(cov_wave=='Second'), pid+100000, pid)) %>%
    arrange(start) %>%
    mutate(n_pid = as.numeric(factor(pid, levels=unique(alog_all_spec$pid)))) %>%
    arrange(n_pid, start)

  ggplot(alog_plot, aes(x=start, y=n_pid, color=cov_wave)) +
    geom_point(size=0.8, alpha=0.5) +
    labs(title='Dotted chart representing the complete COVID-19 rehabilitation log.',
         x='Timestamp',
         y='Patients',
         color='COVID-19 Wave') +
    scale_x_date(date_breaks='3 months', date_labels='%Y-%m') +
    theme_minimal() +
    theme(text=element_text(size=10),
          plot.title=element_text(size=10, face='bold'),
          axis.title=element_text(size=10, face='bold'))
  ggsave(out_path, bg='white', units='in', width=6.25,
         height=4.75, dpi=1200)
}


