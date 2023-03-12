test_that("subset_sig_by_step works", {
  ## test process_data()
  proc_data <- process_data(mastR::im_data_6,
                            group_col = "celltype:ch1",
                            target_group = "NK")
  expect_true(is(proc_data, "DGEList"))

  ## test plot_diagnostics()
  expect_null(plot_diagnostics(proc_data$counts, proc_data$vfit$E,
                               group_col = proc_data$samples$group))

  ## test plot_mean_var()
  plot_mean_var(proc_data)
  expect_null(plot_mean_var(proc_data))

  ## test select_sig()
  sig <- select_sig(proc_data$tfit)
  expect_true(is(sig, "GeneSetCollection"))
})